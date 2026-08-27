"""
    Advection

One-dimensional advection schemes, as types.

Each scheme is an immutable value; `advect!` dispatches on it. Scratch memory,
where a scheme needs any, is an explicit argument obtained from `workspace`.

This replaces the previous `generate_solver` closure factories. The reasons,
in order of weight:

  * **Reentrancy, and so the 2D2P goal.** A 2D2P step sweeps O(N²) independent
    lines per direction, which is work you want to hand to `Threads.@threads`.
    The old design could not: `PFCNonUniform` captured one `f_tmp` shared by
    every line, and `SemiLagrangian` held its buffer in the closure. Moving the
    arrays and the scratch into the argument list makes the scheme value
    immutable and shareable, with only the workspace per task.
  * The duplication was never a code-generation problem. It was an argument
    baked into a closure, and it is fixed by passing the argument.
  * Named types have a method table: docstrings, `methods()`, `precompile`
    targets, and readable stack traces. Anonymous closures have none of those.
  * Invalid options fail at construction. The `Symbol`-valued options are
    singleton types now, so `Godunov(:typo)` is a `MethodError` where it is
    written rather than a call on `nothing` from inside the hot loop.
"""
module Advection

using Interpolations

export AbstractAdvection1D, advect!, workspace,
       Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC, PFCNonUniform,
       PiecewiseConstant, PiecewiseLinear, NoLimiter, VanLeer,
       LinearSpline, QuadraticSpline, CubicSpline

# ---------------------------------------------------------------- option types

"""Reconstruction of the cell interface used by [`Godunov`](@ref)."""
abstract type AbstractReconstruction end

"""Piecewise-constant reconstruction. Reduces `Godunov` to first-order upwind."""
struct PiecewiseConstant <: AbstractReconstruction end

"""
Piecewise-linear reconstruction.

!!! warning
    Unconditionally unstable without a limiter: with [`NoLimiter`](@ref) the
    flux collapses to a centred one, and centred flux with forward Euler
    amplifies every mode. Pair it with [`VanLeer`](@ref).
"""
struct PiecewiseLinear <: AbstractReconstruction end

"""Flux limiter. Callable as `limiter(r)`."""
abstract type AbstractLimiter end

"""No limiting: `r -> 1`."""
struct NoLimiter <: AbstractLimiter end

"""Van Leer limiter [Van Leer, J. Comput. Phys., 14 (4), 361 (1974)]."""
struct VanLeer <: AbstractLimiter end

(::NoLimiter)(r) = 1.0
(::VanLeer)(r) = (r + abs(r))/(1.0 + abs(r))

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
Boundaries are periodic throughout.
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

Finite-volume scheme with the given interface reconstruction and flux limiter.
"""
struct Godunov{R<:AbstractReconstruction, L<:AbstractLimiter} <: AbstractAdvection1D
    reconstruction::R
    limiter::L
end
Godunov(r::AbstractReconstruction) = Godunov(r, NoLimiter())

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

`fmin` and `fmax` are required and bracket the data. By Liouville's theorem `f`
stays within the extrema of the *initial* condition along characteristics, so
the caller supplies the global initial bounds; they cannot be derived per line
inside a multidimensional sweep. A default here would be silently wrong for any
data exceeding it, which is a mistake this package has already made once — the
two former overloads disagreed by a factor of 740 because one defaulted
`fmax = 1`.

`checked` guards against exactly that, at the cost of a `minimum`/`maximum` pass
per call: measured at 13% of the step at N = 10000. It is a type parameter, so
`checked = false` compiles the check away entirely for production runs where the
bounds are known good.
"""
struct PFC{T<:AbstractFloat, Checked} <: AbstractAdvection1D
    fmin::T
    fmax::T
end

function PFC(; fmin, fmax, checked::Bool = true)
    lo, hi = promote(float(fmin), float(fmax))
    return PFC{typeof(lo), checked}(lo, hi)
end

"""
    PFCNonUniform(Δx; fmin, fmax)

Positive Flux Conservative scheme on a static non-uniform grid, `Δx` being the
cell widths. Needs a [`workspace`](@ref).

The limiter coefficient is computed per cell triple. A single global value
would let one refined region tighten the limiter across the whole domain.
"""
struct PFCNonUniform{T<:AbstractFloat} <: AbstractAdvection1D
    Δx::Vector{T}
    ξ::Vector{T}
    fmin::T
    fmax::T
end

"""
    slope_limit(r)

PFC limiter coefficient for a cell triple whose smallest-to-largest spacing
ratio is `r ∈ (0, 1]`. Equals 2 on a locally uniform grid.
"""
slope_limit(r) = (1.0 + r)*(1.0 + 2r)/(3.0 + (r - 1.0/r)^2)

function PFCNonUniform(Δx_::AbstractVector{T}; fmin, fmax) where {T<:AbstractFloat}
    Δx = collect(Δx_)
    n = length(Δx)
    ξ = similar(Δx)
    for i in eachindex(Δx)
        d₋ = Δx[i == 1 ? n : i-1]
        d₊ = Δx[i == n ? 1 : i+1]
        ξ[i] = slope_limit(min(d₋, Δx[i], d₊)/max(d₋, Δx[i], d₊))
    end
    fmn, fmx = promote(float(fmin), float(fmax))
    return PFCNonUniform{T}(Δx, ξ, T(fmn), T(fmx))
end

# ------------------------------------------------------------------ workspace

"""
    workspace(scheme, n)

Scratch memory for `scheme` at problem size `n`, or `nothing` when it needs
none. One workspace per task: that is what makes the schemes safe to run
concurrently over the lines of a multidimensional sweep.
"""
workspace(::AbstractAdvection1D, ::Integer) = nothing

struct SplineWorkspace{T}
    buffer::Vector{T}
end
workspace(::SemiLagrangian{LinearSpline}, n::Integer) = SplineWorkspace(Vector{Float64}(undef, n + 1))
workspace(::SemiLagrangian, n::Integer) = SplineWorkspace(Vector{Float64}(undef, n))

struct PFCWorkspace{T}
    accumulator::Vector{T}
end
workspace(s::PFCNonUniform{T}, n::Integer) where {T} = PFCWorkspace(Vector{T}(undef, n))

"""
    advect!(dest, src, scheme, c[, ws])

Advance `src` one step into `dest` at Courant number `c`.

The four-argument form allocates a workspace when the scheme needs one; pass
one explicitly in any loop that runs more than once.
"""
advect!(dest, src, scheme::AbstractAdvection1D, c) =
    advect!(dest, src, scheme, c, workspace(scheme, length(dest)))

include(joinpath("schemes", "upwind.jl"))
include(joinpath("schemes", "lax_wendroff.jl"))
include(joinpath("schemes", "godunov.jl"))
include(joinpath("schemes", "semi_lagrangian.jl"))
include(joinpath("schemes", "pfc.jl"))
include(joinpath("schemes", "pfc_nonuniform.jl"))

end # module
