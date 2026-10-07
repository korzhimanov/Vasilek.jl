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

using Interpolations
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

"""Flux limiter. Callable as `limiter(r)`."""
abstract type AbstractLimiter end

"""No limiting: `r -> 1`."""
struct NoLimiter <: AbstractLimiter end

"""Van Leer limiter [Van Leer, J. Comput. Phys., 14 (4), 361 (1974)]."""
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
(::VanLeer)(r) = (r + abs(r))/(1 + abs(r))
(::Superbee)(r) = max(zero(r), min(2r, one(r)), min(r, oftype(r, 2)))

"""Interpolating spline used by [`SemiLagrangian`](@ref)."""
abstract type AbstractSpline end

"""Linear B-spline. Equivalent to upwind for `0 < c < 1`."""
struct LinearSpline <: AbstractSpline end

"""Quadratic B-spline with periodic boundaries."""
struct QuadraticSpline <: AbstractSpline end

"""Cubic B-spline with periodic boundaries."""
struct CubicSpline <: AbstractSpline end

_bspline(::LinearSpline) = BSpline(Linear())
_bspline(::QuadraticSpline) = BSpline(Quadratic(Periodic(OnCell())))
_bspline(::CubicSpline) = BSpline(Cubic(Periodic(OnCell())))

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
`|c| ≤ 1`. At `|c| = 1` the correction vanishes and the step is an exact
one-cell shift.

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

With `checked = true` every call verifies `fmin ≤ minimum(src)` and
`maximum(src) ≤ fmax` and throws a `DomainError` otherwise, at the cost of one
`minimum`/`maximum` pass. `checked` is a type parameter, so `checked = false`
compiles the check away for runs whose bounds are known good.
"""
struct PFC{T<:AbstractFloat, Checked} <: AbstractAdvection1D
    fmin::T
    fmax::T
    # `Checked` is the `checked` flag, and `advect!` branches on it: anything
    # but a Bool used to construct, and then failed inside the first step.
    # The bracket is checked here, which every PFC is built through, rather
    # than in the keyword form only: `PFC{Float64, false}(2.0, 0.0)` built.
    function PFC{T,Checked}(fmin, fmax) where {T<:AbstractFloat, Checked}
        Checked isa Bool || _err_checked(PFC, Checked)
        lo, hi = T(fmin), T(fmax)
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
[`PFC`](@ref)'s, with the same check.

`advect!` takes a displacement `α = vΔt` here, not a Courant number, since a
non-uniform grid has none, and it is bounded by the *narrowest* cell:
`|α| ≤ minimum(Δx)`. Every cell gives up its outgoing flux alone, so a wider
step is wrong in the narrowest cell however wide its neighbours are.
[`OnGrid`](@ref) splits a longer step into sub-steps.
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
    function PFCNonUniform{T,Checked}(Δx, ξ, Δxmin, fmin, fmax) where {T<:AbstractFloat, Checked}
        Checked isa Bool || _err_checked(PFCNonUniform, Checked)
        _check_widths(Δx)
        length(ξ) == length(Δx) || throw(DimensionMismatch(
            "PFCNonUniform needs a ξ per cell, got $(length(ξ)) for $(length(Δx)) cells"))
        T(Δxmin) == minimum(Δx) || throw(ArgumentError(
            "PFCNonUniform: Δxmin = $Δxmin is not minimum(Δx) = $(minimum(Δx))"))
        lo, hi = T(fmin), T(fmax)
        _check_bracket(lo, hi)
        return new{T,Checked}(Δx, ξ, Δxmin, lo, hi)
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
    fmn, fmx = promote(float(fmin), float(fmax))
    # Outward, so that a narrower `T` keeps the data the bounds admit: a Float32
    # grid's `fmax = maximum(f)` of Float64 data rounded below it as often as not,
    # and the first checked step refused the data's own peak.
    return PFCNonUniform{T, checked}(Δx, ξ, minimum(Δx), T(fmn, RoundDown), T(fmx, RoundUp))
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
concurrently over the lines of a multidimensional sweep.
"""
workspace(::AbstractAdvection1D, ::Integer, ::Type = Float64) = nothing

struct SplineWorkspace{T}
    buffer::Vector{T}
end
workspace(::SemiLagrangian{LinearSpline}, n::Integer, ::Type{T} = Float64) where {T} =
    SplineWorkspace(Vector{T}(undef, n + 1))
workspace(::SemiLagrangian, n::Integer, ::Type{T} = Float64) where {T} =
    SplineWorkspace(Vector{T}(undef, n))

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
    wrong as a shorter one). A scheme that takes none now refuses any other,
    which it used to accept and ignore; see [`_validate_workspace`](@ref);
  * a workspace whose buffer shares memory with `dest` or `src`
    (`Base.mightalias` again): the spline prefilter overwrites its buffer,
    which `dest` is then sampled from, and `PFCNonUniform` accumulates into
    its own while it still reads `src`;
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

"""
    _validate_workspace(scheme, ws, n)

Check that `ws` is what `workspace(scheme, n)` would have returned. Exact
equality, not a lower bound: the spline prefilter runs over the whole buffer,
so a longer one is as wrong as a shorter one.

A scheme takes `nothing` unless it has a method here for its own workspace, so
a scheme that needs none refuses whatever else it is given: `Upwind` took
another scheme's workspace without a word, which is how a sweep that mixes up
its workspaces goes unnoticed until a scheme that reads one gets the wrong one.
A new scheme with a workspace adds its method, and one to `_scratch_aliases`,
or is refused loudly.
"""
_validate_workspace(::AbstractAdvection1D, ::Nothing, ::Integer) = nothing
_validate_workspace(scheme::AbstractAdvection1D, ws, ::Integer) = _err_foreign_workspace(scheme, ws)
_validate_workspace(::SemiLagrangian{LinearSpline}, ws::SplineWorkspace, n::Integer) =
    length(ws.buffer) == n + 1 ? nothing : _err_workspace(length(ws.buffer), n + 1)
_validate_workspace(::SemiLagrangian, ws::SplineWorkspace, n::Integer) =
    length(ws.buffer) == n ? nothing : _err_workspace(length(ws.buffer), n)
function _validate_workspace(p::PFCNonUniform, ws::PFCWorkspace, n::Integer)
    length(p.Δx) == n || _err_length(n, length(p.Δx))
    length(ws.accumulator) == n || _err_workspace(length(ws.accumulator), n)
    return nothing
end

# Whether the buffer a workspace writes shares memory with `dest` or `src`, for
# `_validate`. A scheme with a workspace adds its method here as well as to
# `_validate_workspace`: a workspace type without one is a MethodError, rather
# than taken as sharing nothing.
_scratch_aliases(::Nothing, dest, src) = false
_scratch_aliases(ws::SplineWorkspace, dest, src) =
    Base.mightalias(ws.buffer, dest) || Base.mightalias(ws.buffer, src)
_scratch_aliases(ws::PFCWorkspace, dest, src) =
    Base.mightalias(ws.accumulator, dest) || Base.mightalias(ws.accumulator, src)

"""
    _check_bounds(src, fmin, fmax)

The check `PFC` and `PFCNonUniform` run on every call when built with
`checked = true`: that `src` lies inside the `[fmin, fmax]` their limiters are
built on, or a `DomainError`. Both call this one function, so they refuse the
same data with the same message. Inlined, so each scheme's `if Checked` removes
it whole. An exception rather than `@assert`, which is a debugging aid that may
be compiled out, and this is part of the schemes' contract.
"""
@inline function _check_bounds(src, lo, hi)
    m = minimum(src)
    lo ≤ m || _err_bound_low(lo, m)
    M = maximum(src)
    M ≤ hi || _err_bound_high(hi, M)
    return nothing
end

@noinline _err_bound_low(lo, m) = throw(DomainError(m,
    "fmin = $lo exceeds minimum(src) = $m: the PFC limiter is built on [fmin, fmax]"))
@noinline _err_bound_high(hi, M) = throw(DomainError(M,
    "fmax = $hi is below maximum(src) = $M: the PFC limiter is built on [fmin, fmax]"))

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
right, and `PFC` loses positivity first. At `|c| = 1` each is an exact one-cell
shift.

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
    OnGrid(scheme, Δz)
    OnGrid(scheme::PFCNonUniform[, Δz])

`scheme` on the grid whose cell widths are `Δz`. An `OnGrid` is itself a scheme,
and its fourth `advect!` argument is always a **displacement**, a length,
whatever `scheme` takes:

  * a scheme that takes a Courant number gets the displacement divided by the
    grid's spacing, and is refused on a non-uniform grid (each width has to
    match the first to `1e-12` of it), rather than given one of its spacings
    and left wrong by the ratio between them -- an `ArgumentError`;
  * `PFCNonUniform` takes the displacement as it is. It carries its grid, so
    `Δz` may be left out; given, it must be that grid.

An `OnGrid` given again with a grid, `OnGrid(og, Δz)`, is `og` if `Δz` is its
grid to the same `1e-12`, and refused otherwise.

A step wider than a cell -- the narrowest one, for `PFCNonUniform` -- is split
into the fewest equal sub-steps that fit, [`substeps`](@ref), since `advect!`
refuses it whole and a translation by `α` is `m` translations by `α/m`. A
scheme with no Courant limit, `SemiLagrangian`, takes it in one step; see
[`nsubsteps`](@ref). Every step a scheme takes on its own is one step here, to
the bit: within a cell `m = 1`, and `α/1` is `α`.

`workspace(og, n[, T])` is scratch for the scheme and a line buffer for the
sub-steps, in the grid's element type unless `T` is given; `n` is the grid's
length.
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

function OnGrid(s::PFCNonUniform, Δz = s.Δx)
    length(Δz) == length(s.Δx) || throw(DimensionMismatch(
        "PFCNonUniform has $(length(s.Δx)) cells, the grid $(length(Δz))"))
    all(isapprox(d, w; rtol = 1e-12) for (d, w) in zip(Δz, s.Δx)) || throw(ArgumentError(
        "PFCNonUniform carries cell widths other than the grid's: it would advect " *
        "on its own grid while the run sums on this one"))
    return OnGrid(s, s.Δxmin, length(s.Δx))
end

function OnGrid(s::AbstractAdvection1D, Δz)
    isempty(Δz) && throw(ArgumentError("OnGrid needs a grid of at least one cell"))
    h = float(first(Δz))
    all(d -> isapprox(d, h; rtol = 1e-12), Δz) || throw(ArgumentError(
        "$(nameof(typeof(s))) takes a Courant number, which a non-uniform grid " *
        "does not have (spacings range over $(extrema(Δz))); use PFCNonUniform here"))
    return OnGrid(s, h, length(Δz))
end

# Already on a grid: kept, if it is this one, to the tolerance the grid is
# taken at. Exactly was too strict: `L/N` and the widths of the cell centres
# `(j - 1/2)L/N` differ in the last bit for most `N`, and refused an `OnGrid`
# built for the very grid it was given with.
function OnGrid(og::OnGrid, Δz)
    same = OnGrid(og.scheme, Δz)
    isapprox(same.h, og.h; rtol = 1e-12) && same.n == og.n || throw(ArgumentError(
        "an OnGrid on $(og.n) cells of $(og.h) given for a grid of $(same.n) cells of $(same.h)"))
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

function workspace(og::OnGrid{S,T}, n::Integer, ::Type{U} = T) where {S, T, U}
    n == og.n || _err_grid(og, n)
    return OnGridWorkspace(workspace(og.scheme, n, U), Vector{U}(undef, n))
end

@noinline _err_grid(og, n) = throw(DimensionMismatch(
    "$(nameof(typeof(og.scheme))) is on a grid of $(og.n) cells, not $n"))

"""
    substeps(α, h)

The fewest equal parts of `α` that are each no longer than `h`, and one for an
`α` that is not finite, which `advect!` then refuses itself.
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
wherever the step alone would have been taken.
"""
nsubsteps(og::OnGrid, α) = substeps(α, og.h)
nsubsteps(::OnGrid{<:SemiLagrangian}, α) = 1

# The scheme's own fourth argument for a sub-step of displacement `a`: a
# Courant number, `a/h` -- divided, not multiplied by `1/h`, which rounds
# differently -- or, for `PFCNonUniform`, `a` itself.
_argument(og::OnGrid, a) = a/og.h
_argument(::OnGrid{<:PFCNonUniform}, a) = a

function advect!(dest, src, og::OnGrid, α, ws::OnGridWorkspace)
    length(src) == og.n || _err_grid(og, length(src))
    length(ws.buffer) == og.n || _err_workspace(length(ws.buffer), og.n)
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
_scratch_aliases(ws::OnGridWorkspace, dest, src) =
    Base.mightalias(ws.buffer, dest) || Base.mightalias(ws.buffer, src)

include(joinpath("schemes", "upwind.jl"))
include(joinpath("schemes", "lax_wendroff.jl"))
include(joinpath("schemes", "godunov.jl"))
include(joinpath("schemes", "semi_lagrangian.jl"))
include(joinpath("schemes", "pfc.jl"))
include(joinpath("schemes", "pfc_nonuniform.jl"))

end # module
