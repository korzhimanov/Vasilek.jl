"""
    Advection

One-dimensional advection schemes, as types.

Each scheme is an immutable value and `advect!(dest, src, scheme, c, ws)`
dispatches on it. Scratch memory, where a scheme needs any, is an explicit
argument from [`workspace`](@ref), one per task, so a scheme value can be
shared by every line of a multidimensional sweep and by every thread. Options
are types (`PiecewiseLinear()`, `VanLeer()`, …), so an invalid one fails where
it is written. Boundaries are periodic.
"""
module Advection

import ..workspace

export AbstractAdvection1D, advect!, workspace,
       Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC, PFCNonUniform, OnGrid,
       PiecewiseConstant, PiecewiseLinear, NoLimiter, VanLeer, Superbee,
       LinearSpline, QuadraticSpline, CubicSpline

# ---------------------------------------------------------------- option types

"""Reconstruction of the cell interface used by [`Godunov`](@ref)."""
abstract type AbstractReconstruction end

"""Piecewise-constant reconstruction. Reduces `Godunov` to first-order upwind."""
struct PiecewiseConstant <: AbstractReconstruction end

"""
Piecewise-linear reconstruction, its slope set by the limiter. With
[`NoLimiter`](@ref) the slope is the downwind difference and `Godunov` reduces
to [`LaxWendroff`](@ref), second order and not monotone; with
[`VanLeer`](@ref) or [`Superbee`](@ref) it is total-variation diminishing up to
`|c| = 1`.
"""
struct PiecewiseLinear <: AbstractReconstruction end

"""
Flux limiter. Callable as `limiter(r)` for any `r`, `±Inf` included:
[`Godunov`](@ref) passes `r = ±Inf` wherever its ratio of differences exceeds
`floatmax`, so a limiter of one's own must return a finite value there.
"""
abstract type AbstractLimiter end

"""No limiting: `r -> 1`."""
struct NoLimiter <: AbstractLimiter end

"""
Van Leer limiter [Van Leer, J. Comput. Phys., 14 (4), 361 (1974)]:
`r -> (r + |r|)/(1 + |r|)`, which is 0 for `r ≤ 0` and tends to 2 as `r` grows;
it is 2 at `r = Inf`.
"""
struct VanLeer <: AbstractLimiter end

"""
Superbee limiter [Roe, Annu. Rev. Fluid Mech. 18, 337 (1986)]:
`r -> max(0, min(2r, 1), min(r, 2))`, the upper edge of Sweby's second-order
total-variation-diminishing region and so the most compressive limiter that
keeps [`Godunov`](@ref) with [`PiecewiseLinear`](@ref) second order and TVD up
to `|c| = 1`.

Use it for sharp fronts. Its compression is anti-diffusion: on smooth data it
steepens slopes, so the L² norm of a perturbation can grow, which biases a
damping rate low in a Vlasov–Poisson run. For smooth phase-space dynamics prefer
[`VanLeer`](@ref) or [`PFC`](@ref). `test/test_comparison.jl` and
`verification/scheme-comparison.jl` quantify both.
"""
struct Superbee <: AbstractLimiter end

(::NoLimiter)(r) = one(r)

# Evaluated as written, `(r + |r|)/(1 + |r|)` overflowed: `r + |r|` is Inf from
# floatmax/2 up, and at r = ±Inf the quotient is Inf/Inf or NaN/Inf, NaN.
# `_ratio` gives such an r wherever |Δ⁻/Δ⁺| exceeds floatmax -- in Float64, say,
# a subnormal downwind difference beside an upwind one of order one, in Float16
# one merely 65504 times smaller -- and the step came out Inf or NaN there.
#
# So r is held to floatmax/2, where the formula is already exactly 2, and a
# negative r, for which it is 0, goes to 0 before any arithmetic. Wherever the
# formula was finite this is it to the bit: `q + q` is the `r + |r|` of a
# positive r, and a negative one still gives +0. NaN stays NaN. `ifelse` for the
# reason `_ratio` gives, in schemes/godunov.jl.
#
# Clamping r in `_ratio` instead would serve any limiter that can take
# floatmax/2, but it measured 6% to 19% over the unguarded VanLeer step, and up
# to 12% on Superbee, which takes ±Inf as it is; this form measured at most 14%
# over it. Julia 1.10 and 1.13, Float32 and Float64, best of 60 interleaved runs
# on a shared host.
function (::VanLeer)(r::AbstractFloat)
    h = floatmax(r)/2
    q = ifelse(r ≤ 0, zero(r), ifelse(r > h, h, r))
    return (q + q)/(1 + q)
end
# A type without `floatmax`, a rational or a dual number, is off the vectorized
# path and takes the same guards with a branch: r ≤ 0 gives 0, and an infinite
# 2r the limit 2, with a zero derivative.
function (::VanLeer)(r)
    q = ifelse(r ≤ 0, zero(r), r)
    isinf(q + q) && return oftype(one(q)/one(q), 2)
    return (q + q)/(1 + q)
end
(::Superbee)(r) = max(zero(r), min(2r, one(r)), min(r, oftype(r, 2)))

"""Interpolating spline used by [`SemiLagrangian`](@ref)."""
abstract type AbstractSpline end

"""Linear B-spline. Equivalent to upwind for `0 < c < 1`."""
struct LinearSpline <: AbstractSpline end

"""Quadratic B-spline with periodic boundaries."""
struct QuadraticSpline <: AbstractSpline end

"""Cubic B-spline with periodic boundaries."""
struct CubicSpline <: AbstractSpline end

# The diagonal and off-diagonal `(a, b)` of a quadratic or cubic B-spline's
# periodic prefilter, in `T`: the spline's weights on a knot and on each of its
# neighbours, evaluated at that knot.
_diagonals(::Type{T}, ::QuadraticSpline) where {T} = (T(3)/T(4), T(1)/T(8))
_diagonals(::Type{T}, ::CubicSpline) where {T} = (T(2)/T(3), T(1)/T(6))

# --------------------------------------------------------------- scheme types

"""
    AbstractAdvection1D

A one-dimensional advection scheme. Advance one step with

    advect!(dest, src, scheme, c)
    advect!(dest, src, scheme, c, ws)

where `c` is the Courant number and `ws` is scratch from [`workspace`](@ref).
Boundaries are periodic throughout. Every scheme but `SemiLagrangian` needs
`|c| ≤ 1`, and `advect!` refuses more; see [`_validate_courant`](@ref).

The fourth argument is a Courant number for every scheme but `PFCNonUniform`,
which takes a displacement, since a non-uniform grid has no single Courant
number. [`OnGrid`](@ref) hides the difference: `OnGrid(scheme, Δz)` is the
scheme on the grid of cell widths `Δz`, takes a displacement whatever `scheme`
takes, refuses a uniform scheme on a non-uniform grid, and splits a step wider
than a cell.
"""
abstract type AbstractAdvection1D end

"""First-order upwind. Total-variation diminishing, first order."""
struct Upwind <: AbstractAdvection1D end

"""
Lax–Wendroff. Second order, and **not** monotone — it overshoots at
discontinuities, as Godunov's theorem requires of any linear scheme above
first order.
"""
struct LaxWendroff <: AbstractAdvection1D end

"""
    Godunov(reconstruction, limiter = NoLimiter())

Finite-volume scheme with the given interface reconstruction and flux limiter:
reconstruct `f` in each cell, carry the reconstruction exactly for one step,
and average it back onto the cells.

The flux through an interface is the upwind cell's reconstruction averaged over
the strip that crosses the interface in one step. For `PiecewiseLinear` and
`c > 0` (the other direction is the mirror image) the strip's midpoint lies
`cΔx/2` upwind of the interface, which is `(1 − c)Δx/2` downwind of the upwind
cell's centre, so

    Φᵢ₋½ = c·(fᵢ₋₁ + φ(r)(1 − c)(fᵢ − fᵢ₋₁)/2),    r = (fᵢ₋₁ − fᵢ₋₂)/(fᵢ − fᵢ₋₁),

which is Sweby's flux-limited Lax–Wendroff [Sweby, SIAM J. Numer. Anal. 21 (5),
995 (1984)]. `φ = 1` gives [`LaxWendroff`](@ref) and `φ = 0` gives
[`Upwind`](@ref). [`VanLeer`](@ref) and [`Superbee`](@ref) keep `φ(r) ≤ 2` and
`φ(r)/r ≤ 2`, which makes the scheme total-variation diminishing for every
`|c| ≤ 1`. At `|c| = 1` the correction vanishes and the step is a one-cell
shift to within an ulp; see [`_validate_courant`](@ref).

A cell takes the difference of the fluxes through its two faces,
`f + (Φin − Φout)`, so a constant line stays constant to the bit.

The `(1 − c)` factor is what makes the scheme second order; 0.1's
`:Riemann_linear` lacked it and was first order, and TVD only to `|c| ≤ 1/2`.
A limiter given with `PiecewiseConstant()` would have no slope to act on and is
refused.
"""
struct Godunov{R<:AbstractReconstruction, L<:AbstractLimiter} <: AbstractAdvection1D
    reconstruction::R
    limiter::L
end
Godunov(r::AbstractReconstruction) = Godunov(r, NoLimiter())

# `PiecewiseConstant` has no slope to limit, so a limiter given with it would be
# silently ignored. Refused instead.
function Godunov(r::PiecewiseConstant, l::AbstractLimiter)
    l isa NoLimiter || throw(ArgumentError(
        "Godunov(PiecewiseConstant(), $(nameof(typeof(l)))()): a piecewise-constant " *
        "reconstruction has no slope, so the limiter would be ignored; use " *
        "PiecewiseLinear() to limit one"))
    return Godunov{PiecewiseConstant, typeof(l)}(r, l)
end

"""
    SemiLagrangian(spline = CubicSpline())

Backward characteristic tracing with B-spline interpolation. Needs a
[`workspace`](@ref). Global order equals the spline degree.
The quadratic and cubic splines prefilter their data, [`_prefilter!`](@ref).
"""
struct SemiLagrangian{S<:AbstractSpline} <: AbstractAdvection1D
    spline::S
end
SemiLagrangian() = SemiLagrangian(CubicSpline())

"""
    PFC(; fmin, fmax, checked = true)

Positive Flux Conservative scheme on a uniform grid.

`fmin` and `fmax` are required and must bracket the data: the limiter is built
on them, and data outside them is corrupted rather than clipped. By Liouville's
theorem `f` stays within the extrema of the *initial* condition, so pass the
global initial bounds; they cannot be derived per line inside a sweep. There is
deliberately no default.

With `checked = true` every call verifies that `src` lies inside
`[fmin, fmax]` to round-off and throws a `DomainError` otherwise, and clamps
what it returns into `[fmin, fmax]`, at the cost of one `minimum`/`maximum`
pass and a clamp. The limiter keeps the data inside in exact arithmetic, and
rounding took a cell an ulp out -- at `|c| = 1`, a one-cell shift only to within
an ulp (see [`_validate_courant`](@ref)), a few ulps short of it, and where a
plateau at a bound met rough data -- which the next checked call refused: under
`fmax = 0.987`, `c = 1` took `[0.195, 0.694, 0.566, 0.987]` to a first cell of
`0.9870000000000001`. The clamp cuts that ulp, and the mass with it; a line
inside is returned as the arithmetic gives it, to the bit. Another scheme's
step leaves such ulps too, and the check takes them: round-off is
[`_roundoff`](@ref) past a bound, and the check is [`_check_bounds`](@ref).
What it takes comes back clamped too, even at `c = 0`, so that round-off
another scheme leaves does not build up through it from call to call.
`checked` is a type parameter, so `checked = false` compiles both away for runs
whose bounds are known good.
"""
struct PFC{T<:AbstractFloat, Checked} <: AbstractAdvection1D
    fmin::T
    fmax::T
    # `Checked` is the `checked` flag, and `advect!` branches on it: anything
    # but a Bool used to construct, and then failed inside the first step.
    # The bracket is checked here, which every PFC is built through, rather
    # than in the keyword form only: `PFC{Float64, false}(2.0, 0.0)` built.
    # Rounded outward into `T`, `_outward`, so that a narrower `T` keeps the data
    # the bounds admit: a Float32 PFC built on ±0.4029317030250849 rounded both
    # inside it, and its first checked step refused the data's own extrema.
    function PFC{T,Checked}(fmin, fmax) where {T<:AbstractFloat, Checked}
        Checked isa Bool || _err_checked(PFC, Checked)
        lo, hi = _outward(T, fmin, fmax)
        _check_bracket(lo, hi)
        return new{T,Checked}(lo, hi)
    end
end

function PFC(; fmin, fmax, checked::Bool = true)
    lo, hi = promote(float(fmin), float(fmax))
    return PFC{typeof(lo), checked}(lo, hi)
end

"""
    PFCNonUniform(Δx; fmin, fmax, checked = true)

Positive Flux Conservative scheme on a static non-uniform grid, `Δx` being the
cell widths. Needs a [`workspace`](@ref).

The limiter coefficient is computed per cell triple, so a refined region does
not tighten the limiter elsewhere. `fmin`, `fmax` and `checked` are
[`PFC`](@ref)'s, with the same check and the same clamp.

`advect!` takes a displacement `α = vΔt` here, not a Courant number, since a
non-uniform grid has none, and it is bounded by the *narrowest* cell:
`|α| ≤ minimum(Δx)`. Every cell gives up its outgoing flux alone, so a wider
step is wrong in the narrowest cell however wide its neighbours are.
[`OnGrid`](@ref) splits a longer step into sub-steps.

At `|α| = minimum(Δx)` the narrowest cells take a step of exactly one cell,
which rounding took out of the bounds as it did `PFC`'s: with `fmin = 0` and
data that touches it, `α = 0.1` over `fill(0.1, 6)` took
`[0, 0.01268125, 0.01268125, 0, 0, 0.01268125]` to `-1.7e-18` in its second
cell. A checked one clamps it to 0.
"""
struct PFCNonUniform{T<:AbstractFloat, Checked} <: AbstractAdvection1D
    Δx::Vector{T}
    ξ::Vector{T}
    Δxmin::T
    fmin::T
    fmax::T
    # As for `PFC`: `Checked` is the `checked` flag, a Bool, and the invariants
    # are checked here, which every PFCNonUniform is built through. A `Δxmin`
    # wider than the narrowest cell let `advect!` take a step that crossed it.
    # The widths are checked in `T`, as the scheme holds them: `Δxmin` rounded
    # to `T` against widths as given refused Float64 widths of 0.05 for a
    # Float32 scheme, and a width can round to zero or overflow in `T`. The
    # bounds round outward, as `PFC`'s do.
    function PFCNonUniform{T,Checked}(Δx, ξ, Δxmin, fmin, fmax) where {T<:AbstractFloat, Checked}
        Checked isa Bool || _err_checked(PFCNonUniform, Checked)
        Δx, ξ = convert(Vector{T}, Δx), convert(Vector{T}, ξ)
        _check_widths(Δx)
        length(ξ) == length(Δx) || throw(DimensionMismatch(
            "PFCNonUniform needs a ξ per cell, got $(length(ξ)) for $(length(Δx)) cells"))
        T(Δxmin) == minimum(Δx) || throw(ArgumentError(
            "PFCNonUniform: Δxmin = $Δxmin is not minimum(Δx) = $(minimum(Δx))"))
        lo, hi = _outward(T, fmin, fmax)
        _check_bracket(lo, hi)
        return new{T,Checked}(Δx, ξ, T(Δxmin), lo, hi)
    end
end

"""
    slope_limit(r)

PFC limiter coefficient for a cell triple whose smallest-to-largest spacing
ratio is `r ∈ (0, 1]`. Equals 2 on a locally uniform grid.
"""
slope_limit(r) = (1 + r)*(1 + 2r)/(3 + (r - 1/r)^2)

# `@constprop :aggressive` so that `checked`, a value, still reaches the type.
# Without it inference did not carry the default `true` through this constructor,
# on 1.10 or 1.13, and `PFCNonUniform(Δx; fmin, fmax)` inferred only as
# `PFCNonUniform{Float64}`, where before the flag it was concrete.
Base.@constprop :aggressive function PFCNonUniform(Δx_::AbstractVector{T}; fmin, fmax,
                                                   checked::Bool = true) where {T<:AbstractFloat}
    Δx = collect(Δx_)
    n = length(Δx)
    _check_widths(Δx)
    ξ = similar(Δx)
    for i in eachindex(Δx)
        d₋ = Δx[i == 1 ? n : i-1]
        d₊ = Δx[i == n ? 1 : i+1]
        ξ[i] = slope_limit(min(d₋, Δx[i], d₊)/max(d₋, Δx[i], d₊))
    end
    # The inner constructor rounds the bounds outward into `T`: a Float32 grid's
    # `fmax = maximum(f)` of Float64 data rounded below it as often as not, and
    # the first checked step refused the data's own peak.
    fmn, fmx = promote(float(fmin), float(fmax))
    return PFCNonUniform{T, checked}(Δx, ξ, minimum(Δx), fmn, fmx)
end

function _check_widths(Δx)
    n = length(Δx)
    n ≥ 3 || throw(ArgumentError("PFCNonUniform needs at least 3 cells, got $n"))
    all(d -> isfinite(d) && d > 0, Δx) || throw(ArgumentError(
        "PFCNonUniform needs finite, positive cell widths; got extrema $(extrema(Δx))"))
    return nothing
end

@noinline _check_bracket(lo, hi) = lo ≤ hi || throw(ArgumentError(
    "fmin = $lo exceeds fmax = $hi"))
@noinline _err_checked(scheme, c) = throw(ArgumentError(
    "$scheme{T, Checked}: Checked is the `checked` flag and must be true or false, got $(repr(c))"))

# ------------------------------------------------------------------ workspace

"""
    workspace(scheme, n[, T])

Scratch memory for `scheme` at problem size `n` and element type `T`
(`Float64` unless the scheme carries its own), or `nothing` when it needs
none. One workspace per task: that is what makes the schemes safe to run
concurrently over the lines of a multidimensional sweep. `advect!` refuses one
narrower than the data, which the step would round to it.
"""
workspace(::AbstractAdvection1D, ::Integer, ::Type = Float64) = nothing

# Scratch for `SemiLagrangian{S}`: the line's spline `coefficients`, rewritten
# every step (`n + 1` of them for `LinearSpline`, the data and the first point
# again), and the prefilter's factorisation `rdiag` and `z`, built with the
# workspace and only read after, empty for `LinearSpline`; see `_prefilter!`.
# They belong to one spline, hence the `S`, the spline's own type whatever the
# scheme's parameter says. `T` is a float, `float(T)` for any other: the step
# computes in it, and the factorisation holds reciprocals.
struct SplineWorkspace{T, S<:AbstractSpline}
    coefficients::Vector{T}
    rdiag::Vector{T}
    z::Vector{T}
end

workspace(s::SemiLagrangian, n::Integer, ::Type{T} = Float64) where {T} =
    _spline_workspace(s.spline, n, float(T))

_spline_workspace(::LinearSpline, n, ::Type{T}) where {T} =
    SplineWorkspace{T, LinearSpline}(Vector{T}(undef, n + 1), T[], T[])

function _spline_workspace(spline::AbstractSpline, n, ::Type{T}) where {T}
    a, b = _diagonals(T, spline)
    rdiag = zeros(T, n)
    z = zeros(T, n)
    # `A + b·uuᵀ` is the periodic matrix only from three nodes up, and advect!
    # refuses fewer; the workspace is still built, for it to refuse.
    if n ≥ 3
        _factorise!(rdiag, a, b)
        z[1] = z[n] = 1
        _thomas!(z, z, rdiag, b)
        z ./= 1 + b*(z[1] + z[n])
    end
    return SplineWorkspace{T, typeof(spline)}(Vector{T}(undef, n), rdiag, z)
end

# The reciprocal pivots of the Thomas factorisation of the tridiagonal `A` with
# `a` on the diagonal, `b` off it and `a − b` in its two corners, as
# `_prefilter!` splits the periodic matrix. Reciprocals, so that a solve
# multiplies where it would divide: both of its sweeps are serial recurrences,
# whose cost is the latency of the chain, and a division costs several times
# what a multiplication does.
function _factorise!(rdiag, a, b)
    n = length(rdiag)
    d = a - b
    rdiag[1] = 1/d
    for i = 2:n
        d = (i == n ? a - b : a) - b*(b*rdiag[i-1])
        rdiag[i] = 1/d
    end
    return rdiag
end

# Solve `A y = x` into `y`, `A` the tridiagonal whose reciprocal pivots
# `_factorise!` put in `rdiag`; `x` may be `y` itself. The forward sweep reads
# `x` and writes `y`, so the prefilter needs no copy of the data first. Needs
# `x`, `y` and `rdiag` one-based and of one length, which `_validate` and
# `_validate_workspace` check before every step.
function _thomas!(y, x, rdiag, b)
    n = length(y)
    n ≥ 1 || return y
    @inbounds begin
        yᵢ = convert(eltype(y), x[1])
        y[1] = yᵢ
        for i = 2:n
            yᵢ = x[i] - (b*rdiag[i-1])*yᵢ
            y[i] = yᵢ
        end
        yᵢ *= rdiag[n]
        y[n] = yᵢ
        # `y[i]*rdiag[i] - (b*rdiag[i])*yᵢ₊₁` rather than `(y[i] - b*yᵢ₊₁)*rdiag[i]`:
        # a multiplication and a subtraction on the chain, not three operations
        for i = n-1:-1:1
            yᵢ = y[i]*rdiag[i] - (b*rdiag[i])*yᵢ
            y[i] = yᵢ
        end
    end
    return y
end

struct PFCWorkspace{T}
    accumulator::Vector{T}
end
workspace(s::PFCNonUniform{T}, n::Integer, ::Type{S} = T) where {T, S} =
    PFCWorkspace(Vector{S}(undef, n))

"""
    advect!(dest, src, scheme, c[, ws])

Advance `src` one step into `dest` at Courant number `c`.

The four-argument form allocates a workspace of `src`'s element type when the
scheme needs one; pass one explicitly in any loop that runs more than once.

`dest` and `src` must be distinct arrays of equal length, at least three
elements long; `ws` must be the workspace `workspace(scheme, length(src))`
returns; and `|c| ≤ 1` unless the scheme is `SemiLagrangian` -- for
`PFCNonUniform`, whose `c` is a displacement, `|c| ≤ minimum(Δx)`. Each of
those is checked; see [`_validate`](@ref) for why.
"""
advect!(dest, src, scheme::AbstractAdvection1D, c) =
    advect!(dest, src, scheme, c, workspace(scheme, length(dest), float(eltype(src))))

# ----------------------------------------------------------------- validation

"""
    _validate(dest, src, scheme, c, ws)

Argument check run at the top of every `advect!` method. Each of these used to
give a plausible wrong answer instead of an error:

  * `dest` sharing memory with `src` (`Base.mightalias`, so a view counts):
    the stencils read neighbours of `src` that `dest` has already overwritten;
  * unequal lengths, fewer than three cells, or non-one-based indexing;
  * a workspace not exactly the one `workspace(scheme, length(src))` returns
    (the spline prefilter runs over the whole buffer, so a longer one is as
    wrong as a shorter one, and a workspace built for another spline holds
    another spline's factorisation). A scheme that takes none now refuses any
    other, which it used to accept and ignore; see [`_validate_workspace`](@ref);
  * a workspace narrower than the data, `_validate_eltype`: a step works the
    line in its workspace's type and rounded the data to it, 5.9e-8 off for
    Float64 data under a Float32 grid's `PFCNonUniform`;
  * a workspace whose buffers share memory with `dest` or `src`
    (`Base.mightalias` again): the spline prefilter overwrites its
    coefficients, which `dest` is then sampled from, and `PFCNonUniform`
    accumulates into its own while it still reads `src`;
  * a step the scheme cannot take; see [`_validate_courant`](@ref).

The checks are comparisons outside the loop and allocate nothing; the error
paths are `@noinline`.
"""
@inline function _validate(dest, src, scheme::AbstractAdvection1D, c, ws)
    Base.require_one_based_indexing(dest, src)
    Base.mightalias(dest, src) && _err_alias()
    length(dest) == length(src) || _err_length(length(dest), length(src))
    length(src) ≥ 3 || _err_short(length(src))
    _validate_workspace(scheme, ws, length(src))
    _validate_eltype(ws, src)
    _scratch_aliases(ws, dest, src) && _err_alias(ws)
    _validate_courant(scheme, c)
    return nothing
end

@noinline _err_alias() = throw(ArgumentError(
    "advect! requires dest and src not to share memory: every scheme reads neighbours of src that " *
    "an aliased dest would already have overwritten"))
@noinline _err_alias(ws) = throw(ArgumentError(
    "advect! requires the workspace not to share memory with dest or src: the scheme overwrites " *
    "its $(nameof(typeof(ws))) buffer while it reads src and writes dest"))
@noinline _err_foreign_workspace(scheme, ws) = throw(ArgumentError(
    "$(nameof(typeof(scheme))) takes no workspace, but was given a $(nameof(typeof(ws))); " *
    "pass workspace(scheme, length(src)), which is nothing for it"))
@noinline _err_length(nd, ns) = throw(DimensionMismatch(
    "advect! requires length(dest) == length(src), got $nd and $ns"))
@noinline _err_short(n) = throw(ArgumentError(
    "advect! needs at least 3 cells, got $n: the schemes' unrolled boundary " *
    "stencils reach two neighbours either side"))
@noinline _err_workspace(got, want) = throw(DimensionMismatch(
    "workspace has length $got but this call needs $want; build it with " *
    "workspace(scheme, length(src))"))
@noinline _err_narrow(T, S) = throw(ArgumentError(
    "the workspace is $T, narrower than the $S data: the step works in the " *
    "workspace's type and would round the data to it; build it with " *
    "workspace(scheme, length(src), $S)"))
@noinline _err_knots(T, n) = throw(ArgumentError(
    "a SemiLagrangian step on $n cells computes its points in $T, which does not " *
    "hold every knot up to $(n + 1) exactly; build the workspace in a wider type, " *
    "workspace(scheme, $n, $(T === Float16 ? Float32 : Float64))"))
@noinline _err_courant_type(T, c) = throw(DomainError(c,
    "the Courant number overflows $T, the type SemiLagrangian computes in"))
@noinline _err_spline(want, S) = throw(ArgumentError(
    "the workspace was built for SemiLagrangian($(nameof(S))()), but this call is " *
    "SemiLagrangian($(nameof(typeof(want)))()); build it with workspace(scheme, length(src))"))

"""
    _validate_workspace(scheme, ws, n)

Check that `ws` is what `workspace(scheme, n)` would have returned. Exact
equality, not a lower bound: the spline prefilter runs over the whole buffer,
so a longer one is as wrong as a shorter one. A spline's workspace also has to
be its own spline's, since it holds that spline's factorisation.

A scheme takes `nothing` unless it has a method here for its own workspace, so
a scheme that needs none refuses whatever else it is given: `Upwind` took
another scheme's workspace without a word, which is how a sweep that mixes up
its workspaces goes unnoticed until a scheme that reads one gets the wrong one.
A new scheme with a workspace adds its method, and one to `_scratch_aliases`,
or is refused loudly.
"""
_validate_workspace(::AbstractAdvection1D, ::Nothing, ::Integer) = nothing
_validate_workspace(scheme::AbstractAdvection1D, ws, ::Integer) = _err_foreign_workspace(scheme, ws)
function _validate_workspace(s::SemiLagrangian, ws::SplineWorkspace, n::Integer)
    # Another spline's workspace has the length and holds a factorisation, of
    # the wrong matrix: the answer it gives looks right and is not. A comparison
    # of types, decided when the method compiles. (A method for the matching
    # pair instead, `SemiLagrangian{S}` with `SplineWorkspace{<:Any, S}`, is not
    # more specific than this one to Julia's dispatch.)
    _spline_of(ws) === typeof(s.spline) || _err_spline(s.spline, _spline_of(ws))
    # `_interior` and `_periodic` need `i - c` and the period exact at every
    # knot: Float16 holds the integers only up to 2048. The comparison is in
    # `T`, which holds `maxintfloat(T)`.
    n + 1 ≤ maxintfloat(_eltype(ws)) || _err_knots(_eltype(ws), n)
    if s.spline isa LinearSpline
        length(ws.coefficients) == n + 1 || _err_workspace(length(ws.coefficients), n + 1)
    else
        length(ws.coefficients) == n || _err_workspace(length(ws.coefficients), n)
        # the prefilter's sweeps index these at the coefficients' length, unchecked
        length(ws.rdiag) == n || _err_workspace(length(ws.rdiag), n)
        length(ws.z) == n || _err_workspace(length(ws.z), n)
    end
    return nothing
end
_spline_of(::SplineWorkspace{<:Any, S}) where {S} = S
_eltype(::SplineWorkspace{T}) where {T} = T
function _validate_workspace(p::PFCNonUniform, ws::PFCWorkspace, n::Integer)
    length(p.Δx) == n || _err_grid(p, length(p.Δx), n)
    length(ws.accumulator) == n || _err_workspace(length(ws.accumulator), n)
    return nothing
end

# A workspace narrower than the data would round it on the way in: the spline
# step computes in its workspace's type, `PFCNonUniform` accumulates in it, and
# `OnGrid`'s sub-steps read the line back from its buffer. A wider one is a
# choice, and allowed. A comparison of types, decided when the method compiles.
_validate_eltype(ws, src) = nothing
_validate_eltype(ws::Union{SplineWorkspace, PFCWorkspace}, src) = _check_narrow(_eltype(ws), src)
_check_narrow(::Type{T}, src) where {T} =
    promote_type(T, float(eltype(src))) === T ? nothing : _err_narrow(T, float(eltype(src)))
_eltype(::PFCWorkspace{T}) where {T} = T

# `c` in the type a spline step computes in. `_validate_courant` has refused a
# non-finite one; one finite in its own type can still overflow a narrower `T`.
function _courant(::Type{T}, c) where {T}
    cT = convert(T, c)
    isfinite(cT) || _err_courant_type(T, c)
    return cT
end

# Whether the buffer a workspace writes shares memory with `dest` or `src`, for
# `_validate`. A scheme with a workspace adds its method here as well as to
# `_validate_workspace`: a workspace type without one is a MethodError, rather
# than taken as sharing nothing.
_scratch_aliases(::Nothing, dest, src) = false
# All three of a spline's vectors: its factorisation is only read during a step,
# but a `dest` on top of it would leave every later step wrong. The linear
# spline has no factorisation, and its two empty vectors nothing to share.
_scratch_aliases(ws::SplineWorkspace, dest, src) =
    _shares(ws.coefficients, dest, src) || _shares(ws.rdiag, dest, src) || _shares(ws.z, dest, src)
_scratch_aliases(ws::SplineWorkspace{<:Any, LinearSpline}, dest, src) =
    _shares(ws.coefficients, dest, src)
_scratch_aliases(ws::PFCWorkspace, dest, src) = _shares(ws.accumulator, dest, src)
_shares(a, dest, src) = Base.mightalias(a, dest) || Base.mightalias(a, src)

"""
    _check_bounds(src, fmin, fmax)

The check `PFC` and `PFCNonUniform` run on every call when built with
`checked = true`: that `src` lies inside the `[fmin, fmax]` their limiters are
built on to round-off, or a `DomainError`. Both call this one function, so they
refuse the same data with the same message. Inlined, so each scheme's
`if Checked` removes it whole. An exception rather than `@assert`, which is a
debugging aid that may be compiled out, and this is part of the schemes'
contract. What a checked call returns it clamps into the bounds, `_hold!`.

Round-off is [`_roundoff`](@ref) past a bound, and only a line outside the
bounds pays for it: one inside is the two comparisons it always was. Another
scheme's step rounds a line past a bound at its own extremum, as one taken
beside the driver's default does, and refusing that ulp stopped the run.
"""
@inline function _check_bounds(src, lo, hi)
    m, M = minimum(src), maximum(src)
    lo ≤ m && M ≤ hi && return nothing
    return _check_roundoff(src, lo, hi, m, M)
end

# A line outside its bounds, cell by cell: each cell past a bound is taken if it
# lies within `_roundoff` of it, at the scale of the bound, the cell and two
# cells either side, as far as one step's stencil reads. Rare, and kept out of
# the inlined check. The refusal names the cell that failed: the line's
# extremum can be one the check took, beside a larger neighbour.
@noinline function _check_roundoff(src, lo, hi, m, M)
    isnan(m) && _err_nan(findfirst(isnan, src))
    F = float(eltype(src))
    i₁, n = firstindex(src), length(src)
    for i in eachindex(src)
        x = src[i]
        x < lo || x > hi || continue
        bound = x < lo ? lo : hi
        S = promote_type(F, float(typeof(bound)))
        scale = isfinite(bound) ? abs(S(bound)) : zero(S)
        for d in -2:2    # an infinite neighbour, under fmax = Inf, is no scale
            y = abs(S(src[i₁ + mod(i - i₁ + d, n)]))
            isfinite(y) && (scale = max(scale, y))
        end
        abs(x - bound) ≤ _roundoff(F, scale) ||
            (x < lo ? _err_bound_low(lo, x, i) : _err_bound_high(hi, x, i))
    end
    return nothing
end

"""
    _roundoff(F, scale)

How far past a bound a value of a line of type `F` is taken for round-off by a
checked `PFC` or `PFCNonUniform`, where `scale` is the largest magnitude among
it, its two neighbours either side and the bound: the values a scheme's stencil
rounded to get there. 1024 ulps of `F` at that scale, 2.3e-13 of it for
Float64 and as many ulps of BigFloat. Nothing that is not finite is round-off.

A line of a type coarser than Float64 is taken only inside its bounds, as
before: no count of its ulps told round-off from data outside them. `Godunov`
with a limiter left a Float32 line 9 ulps over its maximum and a Float16 one
18 within 64 `OnGrid` sub-steps, and 4 ulps took a Float16 line 0.3% below 0.
The driver works its lines in Float64 at least, and never hands one over.

The scale is local, the bound being tested and the cell's neighbours, so that
neither the far bound nor a peak elsewhere on the line widens it: a loose
`fmax = 1e300` took `-5` for round-off past `fmin = 0` while it counted, and a
line's own peak would take a real undershoot of a cell beside zeros. The scale
is taken in a float type both the line and the bound fit, so that a Float32 line
under `fmax = 1e300` does not overflow it and an integer line is judged in
Float64.

A step leaves a Float64 line at most an ulp of that scale past them -- `Upwind`,
the constant `Godunov` and the linear `SemiLagrangian` -- or about an ulp per
sub-step of `OnGrid`, `Godunov` with a limiter, where a plateau at the bound
meets values just below it: 111 ulps over 256 sub-steps and 197 over 1024 at
worst, and none over 4096 or more, where the line has smoothed below the bound.
None of these went below 0 at a jump to it.
"""
function _roundoff(::Type{F}, scale) where {F}
    eps(F) ≤ eps(Float64) && isfinite(scale) || return zero(scale)
    return 1024*eps(F)*scale
end

"""
    _hold!(dest, src, l, u)

`src` clamped into `[l, u]`, into `dest`, which may be `src` itself. `l` and
`u` are the bounds already rounded inward into `dest`'s type, by `_inward`,
so that a narrower `dest` holds nothing past them; a checked `PFC` and
`PFCNonUniform` round them before they write `dest` and clamp what they return
with this. Comparisons and `ifelse`, so that it vectorizes, where `clamp!` cost
a third of a `PFC` step at N = 10000.
"""
function _hold!(dest, src, l, u)
    @inbounds @simd for i in eachindex(dest, src)
        x = src[i]
        dest[i] = ifelse(x > u, u, ifelse(x < l, l, x))
    end
    return dest
end

# A bracket rounded into a type: outward, `_outward`, where it has to hold what
# the wider one held -- `PFCNonUniform`'s bounds, and the driver's `keeps_bounds`
# -- and inward, `_inward`, where nothing outside it may come out of the
# rounding: the clamp. Inward, an `ArgumentError` where no value of `S` lies
# between them: `fmin = fmax = 0.1` rounded to a Float32 `l` above `u`, and the
# clamp put every cell at `u`, below `fmin`. A checked `PFC` and
# `PFCNonUniform` ask before they write.
_outward(::Type{S}, lo, hi) where {S<:AbstractFloat} = (S(lo, RoundDown), S(hi, RoundUp))
function _inward(::Type{S}, lo, hi) where {S<:AbstractFloat}
    l, u = S(lo, RoundUp), S(hi, RoundDown)
    l ≤ u || _err_no_room(S, lo, hi)
    return l, u
end
_inward(::Type{S}, lo, hi) where {S} = (convert(S, lo), convert(S, hi))
@noinline _err_no_room(S, lo, hi) = throw(ArgumentError(
    "no $S lies in [fmin, fmax] = [$lo, $hi], and a checked call returns its line inside " *
    "them: widen the bounds to hold a value of the data's type"))

@noinline _err_bound_low(lo, x, i) = throw(DomainError(x,
    "fmin = $lo exceeds src[$i] = $x by more than round-off: the PFC limiter is " *
    "built on [fmin, fmax]"))
@noinline _err_bound_high(hi, x, i) = throw(DomainError(x,
    "fmax = $hi is below src[$i] = $x by more than round-off: the PFC limiter is " *
    "built on [fmin, fmax]"))
@noinline _err_nan(i) = throw(DomainError(NaN,
    "src[$i] is NaN, which no [fmin, fmax] holds: the PFC limiter is built on them"))

"""
    _validate_courant(scheme, c)

Refuse a step `scheme` cannot take: `|c| > 1` for the schemes with a Courant
limit, and for `PFCNonUniform`, whose `c` is a displacement, `|c|` wider than
its narrowest cell. `c = ±1` exactly is accepted; so is any finite `c` by
`SemiLagrangian`, since characteristic tracing has no Courant limit. `NaN` and
`±Inf` are refused by every scheme.

`Upwind`, `LaxWendroff`, `Godunov` and `PFC` are explicit, and none survives a
characteristic crossing more than one cell per step: past `|c| = 1` round-off
at the grid scale grows every step, slowly enough for the first steps to look
right, and `PFC` loses positivity first. At `|c| = 1` each is a one-cell shift
to within an ulp of the larger of the two cells involved: every flux is then
its upwind cell's whole value, and a cell takes `f[i] + (f[i−1] − f[i])` --
`Upwind`'s `f[i] − (f[i] − f[i−1])` is the same rounding -- which is `f[i−1]`
exactly only where the difference is exact, as it is between neighbours within
a factor of two of each other. `LaxWendroff` rounds more terms. That ulp can
take a cell outside a `PFC`'s bounds, and a checked `PFC` or `PFCNonUniform`
clamps it back.

The bounded method is the default, so a new scheme is checked unless it opts
out, as `SemiLagrangian` does: one that forgets to opt out is refused loudly,
where one that forgot to opt in would be wrong quietly.
"""
_validate_courant(scheme::AbstractAdvection1D, c) =
    abs(c) ≤ 1 ? nothing : _err_courant(scheme, c)
_validate_courant(::SemiLagrangian, c) = isfinite(c) ? nothing : _err_nonfinite(c)
_validate_courant(p::PFCNonUniform, α) =
    abs(α) ≤ p.Δxmin ? nothing : _err_displacement(α, p.Δxmin)

@noinline _err_courant(scheme, c) = throw(DomainError(c,
    "$(nameof(typeof(scheme))) is stable only for |c| ≤ 1: a characteristic may " *
    "not cross more than one cell per step. Take more, smaller steps, or use " *
    "SemiLagrangian, which has no Courant limit"))
@noinline _err_nonfinite(c) = throw(DomainError(c,
    "SemiLagrangian needs a finite Courant number"))
@noinline _err_displacement(α, h) = throw(DomainError(α,
    "PFCNonUniform needs |α| ≤ minimum(Δx) = $h: its fourth argument is a " *
    "displacement, and each cell gives up its outgoing flux alone, so none may " *
    "be crossed in one step. Take more, smaller steps"))

# ------------------------------------------------------------- on a grid

"""
    OnGrid(scheme, Δz; scale = maximum(abs, Δz))
    OnGrid(scheme::PFCNonUniform[, Δz])

`scheme` on the grid whose cell widths are `Δz`. An `OnGrid` is itself a scheme,
and its fourth `advect!` argument is always a **displacement**, a length,
whatever `scheme` takes:

  * a scheme that takes a Courant number gets the displacement divided by the
    grid's spacing, its first width, and is refused on a non-uniform grid,
    rather than given one of its spacings and left wrong by the ratio between
    them -- an `ArgumentError`. Uniform is to rounding, as `make_poisson` takes
    it, [`_nonuniform_at`](@ref): each width matches the first to `1e-10` of
    it, or to 8 ulps of `scale` in `scale`'s type. By default that is the
    widths' own rounding, `scale` the widest of them. Widths taken from
    coordinates carry the coordinates' rounding instead, so pass the largest
    coordinate, in its own type, as `vlasov_poisson` does:
    `OnGrid(scheme, cell_widths(x); scale = maximum(abs, x))`. The grid's
    extent as the default took a Float32 grid of 4096 cells to 2.5e-3 of its
    spacing, and a cell 5e-4 wider than the rest with it;
  * `PFCNonUniform` takes the displacement as it is. It carries its grid, so
    `Δz` may be left out; given, it must be that grid, to `1e-12` of each width
    or 4 ulps of the scheme's own type, so that a Float32 scheme built from a
    Float64 grid's widths takes that grid. `scale` is unused.

An `OnGrid` given again with a grid, `OnGrid(og, Δz)`, is `og` if `Δz` is its
grid to the same rounding, and refused otherwise.

A step wider than a cell -- the narrowest one, for `PFCNonUniform` -- is split
into the fewest equal sub-steps that fit, [`substeps`](@ref), since `advect!`
refuses it whole and a translation by `α` is `m` translations by `α/m`. A
scheme with no Courant limit, `SemiLagrangian`, takes it in one step; see
[`nsubsteps`](@ref). Every step a scheme takes on its own is one step here, to
the bit: within a cell `m = 1`, and `α/1` is `α`.

`workspace(og, n[, T])` is scratch for the scheme and a line buffer for the
sub-steps; `n` is the grid's length. **Its type is the grid's unless `T` is
given, not the data's**: Float64 data on a Float32 grid needs
`workspace(og, n, Float64)`. `advect!` refuses a workspace narrower than the
data, whose sub-steps it rounded to it, 8.9e-8 off at three sub-steps; the
four-argument `advect!` builds one in the data's type.
"""
struct OnGrid{S<:AbstractAdvection1D, T<:AbstractFloat} <: AbstractAdvection1D
    scheme::S
    h::T
    n::Int
    # Every OnGrid is built through here. `h` is what a displacement is divided
    # by, or for `PFCNonUniform` what bounds a sub-step, which has to be its
    # narrowest cell: anything wider let a sub-step cross it.
    function OnGrid{S,T}(scheme, h, n) where {S<:AbstractAdvection1D, T<:AbstractFloat}
        S <: OnGrid && throw(ArgumentError(
            "OnGrid of an OnGrid: the inner one already takes a displacement"))
        isfinite(h) && h > 0 || throw(ArgumentError(
            "OnGrid needs a finite, positive cell width, got h = $h"))
        if scheme isa PFCNonUniform
            n == length(scheme.Δx) && h == scheme.Δxmin || throw(ArgumentError(
                "OnGrid: a PFCNonUniform of $(length(scheme.Δx)) cells, the narrowest " *
                "$(scheme.Δxmin), is on a grid of $n cells with h = $h"))
        end
        return new{S,T}(scheme, h, n)
    end
end
OnGrid(scheme::S, h::T, n::Integer) where {S<:AbstractAdvection1D, T<:AbstractFloat} =
    OnGrid{S,T}(scheme, h, n)

"""
    _nonuniform_at(Δz, scale)

The index of the first spacing in `Δz` that differs from the first by more than
rounding, or `nothing` for a uniform grid: `1e-10` of the first spacing, or 8
ulps of `scale` in `scale`'s own type, whichever is more. A spacing is the
difference of two coordinates and carries their rounding, which scales with
the coordinates, `scale` the largest of them, rather than with the spacing:
Float32 coordinates up to 25 leave widths of π/8 4.9e-6 of themselves apart,
Float64 ones near 1e6 widths of 0.1 1.2e-9 apart, and `OnGrid` refused both
while it took `1e-12` of the spacing. `make_poisson` checks its nodes by the
same rule, [`_within_rounding`](@ref), and [`OnGrid`](@ref) a uniform scheme's
widths with this. NaN is never within rounding.
"""
function _nonuniform_at(Δz, scale)
    h = first(Δz)
    return findfirst(d -> !_within_rounding(d, h, scale), Δz)
end

"""
    _within_rounding(d, h, scale)

Whether the spacing `d` is `h` to rounding, as [`_nonuniform_at`](@ref) takes it.
"""
_within_rounding(d, h, scale) = d == h || abs(d - h) ≤ max(1e-10*abs(h), 8eps(float(scale)))

# `OnGrid`'s default `scale`: the widest width, so that widths given as they are
# are taken to their own rounding. `init`, so that an empty grid reaches
# `OnGrid`'s own refusal.
_widest(Δz) = maximum(abs, Δz; init = zero(float(eltype(Δz))))

# `scale` serves the uniformity check, which a PFCNonUniform does not need: it
# is taken, so that a caller can pass it whatever the scheme. Its widths are
# compared with the grid's in its own type, as `Yee1D` compares its steps: a
# Float32 scheme built from a Float64 grid's widths holds them rounded, and
# `1e-12` refused it on that grid.
function OnGrid(s::PFCNonUniform, Δz = s.Δx; scale = nothing)
    length(Δz) == length(s.Δx) || throw(DimensionMismatch(
        "PFCNonUniform has $(length(s.Δx)) cells, the grid $(length(Δz))"))
    rtol = max(1e-12, 4eps(eltype(s.Δx)))
    all(isapprox(d, w; rtol) for (d, w) in zip(Δz, s.Δx)) || throw(ArgumentError(
        "PFCNonUniform carries cell widths other than the grid's: it would advect " *
        "on its own grid while the run sums on this one"))
    return OnGrid(s, s.Δxmin, length(s.Δx))
end

function OnGrid(s::AbstractAdvection1D, Δz; scale = _widest(Δz))
    isempty(Δz) && throw(ArgumentError("OnGrid needs a grid of at least one cell"))
    h = float(first(Δz))
    _nonuniform_at(Δz, scale) === nothing || throw(ArgumentError(
        "$(nameof(typeof(s))) takes a Courant number, which a non-uniform grid " *
        "does not have (spacings range over $(extrema(Δz))); use PFCNonUniform here"))
    return OnGrid(s, h, length(Δz))
end

# Already on a grid: kept, if it is this one, to the rounding the grid is taken
# to. Exactly was too strict: `L/N` and the widths of the cell centres
# `(j - 1/2)L/N` differ in the last bit for most `N`, and refused an `OnGrid`
# built for the very grid it was given with.
function OnGrid(og::OnGrid, Δz; scale = _widest(Δz))
    same = OnGrid(og.scheme, Δz; scale)
    _within_rounding(same.h, og.h, scale) && same.n == og.n ||
        throw(ArgumentError("an OnGrid on $(og.n) cells of $(og.h) given for a grid of " *
                            "$(same.n) cells of $(same.h)"))
    return og
end

"""
    OnGridWorkspace

Scratch for an [`OnGrid`](@ref): the scheme's own workspace, and a line the
sub-steps after the first read from.
"""
struct OnGridWorkspace{W, T}
    inner::W
    buffer::Vector{T}
end

# The sub-step line holds the data between sub-steps, so it may not be narrower
# than they are; the scheme's own `advect!` checks its own scratch.
_validate_eltype(ws::OnGridWorkspace, src) = _check_narrow(eltype(ws.buffer), src)

function workspace(og::OnGrid{S,T}, n::Integer, ::Type{U} = T) where {S, T, U}
    n == og.n || _err_grid(og.scheme, og.n, n)
    return OnGridWorkspace(workspace(og.scheme, n, U), Vector{U}(undef, n))
end

# A scheme on a grid of `cells` cells, given a line or a workspace of `n`: an
# `OnGrid`'s scheme, or a `PFCNonUniform`, which carries its grid.
@noinline _err_grid(scheme, cells, n) = throw(DimensionMismatch(
    "$(nameof(typeof(scheme))) is on a grid of $cells cells, not $n"))

"""
    substeps(α, h)

The fewest equal parts of `α` that are each no longer than `h`, and one for an
`α` that is not finite, which `advect!` then refuses itself.

A part of exactly `h` is a step of `|c| = 1`, which a checked `PFC` or
`PFCNonUniform` clamps into its bounds: `OnGrid(PFCNonUniform(fill(0.1, 6);
fmin = 0, fmax = 1))` took `α = 0.2` in two steps of exactly 0.1, and the
second refused the `-1.7e-18` the first had left in a cell that held 0, until
the scheme clamped what it returned. For a while those two were given parts
shorter than `h` instead, which a step a few ulps short of a cell still beat.
"""
function substeps(α, h)
    isfinite(α) || return 1
    m = max(1, ceil(Int, abs(α)/h))
    return abs(α/m) > h ? m + 1 : m    # the quotient can round below the integer
end

"""
    nsubsteps(og::OnGrid, α)

How many equal sub-steps [`OnGrid`](@ref) takes the displacement `α` in:
[`substeps`](@ref)`(α, og.h)`, so that none crosses a cell, and 1 for a scheme
with no Courant limit. `SemiLagrangian` has a method returning 1; a scheme of
another type that has none adds one.

`|α/m| ≤ h` is what `PFCNonUniform` checks, and for a scheme that takes a
Courant number it gives `|(α/m)/h| ≤ 1`, division being monotone. Where `α`
fits a cell, `|α| ≤ h` exactly when `|α/h| ≤ 1` after rounding, so `m = 1`
wherever the step alone would have been taken: `Upwind` at `|α| = h` is its own
step at `c = ±1`.
"""
nsubsteps(og::OnGrid, α) = substeps(α, og.h)
nsubsteps(::OnGrid{<:SemiLagrangian}, α) = 1

# The scheme's own fourth argument for a sub-step of displacement `a`: a
# Courant number, `a/h` -- divided, not multiplied by `1/h`, which rounds
# differently -- or, for `PFCNonUniform`, `a` itself.
_argument(og::OnGrid, a) = a/og.h
_argument(::OnGrid{<:PFCNonUniform}, a) = a

function advect!(dest, src, og::OnGrid, α, ws::OnGridWorkspace)
    length(src) == og.n || _err_grid(og.scheme, og.n, length(src))
    length(ws.buffer) == og.n || _err_workspace(length(ws.buffer), og.n)
    _validate_eltype(ws, src)    # whether or not there are sub-steps, for every call alike
    _scratch_aliases(ws, dest, src) && _err_alias(ws)
    m = nsubsteps(og, α)
    a = _argument(og, α/m)
    advect!(dest, src, og.scheme, a, ws.inner)
    for _ in 2:m
        copyto!(ws.buffer, dest)
        advect!(dest, ws.buffer, og.scheme, a, ws.inner)
    end
    return dest
end

# The sub-step line only; the scheme's own `advect!` checks its own scratch.
_scratch_aliases(ws::OnGridWorkspace, dest, src) = _shares(ws.buffer, dest, src)

include(joinpath("schemes", "upwind.jl"))
include(joinpath("schemes", "lax_wendroff.jl"))
include(joinpath("schemes", "godunov.jl"))
include(joinpath("schemes", "semi_lagrangian.jl"))
include(joinpath("schemes", "pfc.jl"))
include(joinpath("schemes", "pfc_nonuniform.jl"))

end # module
