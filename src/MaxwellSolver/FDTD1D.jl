"""
    FDTD1D

The one-dimensional Yee scheme for the transverse fields, as a type.

[`Yee1D`](@ref) is an immutable value holding the step's constants, the
absorbing layers and the source, and `advance!(mesh, op, t, j)` takes one
leapfrog step on a [`YeeMesh1D`](@ref). The operator holds no fields and needs
no scratch, so `workspace(op)` is `nothing` and one value can step any number
of meshes, from any number of tasks.
"""
module FDTD1D

import ..workspace

export YeeMesh1D, PML, Yee1D, advance!

"""
    YeeMesh1D{T}(N)

Staggered fields on `N` cells: `N+1` electric nodes with `N` magnetic cells
between them, all zeroed. `N` is stored as an `Int`.

The two end nodes, `ey[1]`/`ez[1]` and `ey[end]`/`ez[end]`, are **perfect
electric conductor boundaries**. No update touches them — the interior loop
runs `pml.N+2 : Nx-pml.N` and the two absorbing-layer loops stop short of both
ends — so they hold the zero they are constructed with, for all time. They are
not dead storage: the `hz` update reads `ey[end]`, which is how the boundary
condition enters the solution.

Writing a nonzero value into either end node therefore does not seed a wave,
it silently changes the boundary condition to `E = const`. Seed the interior
instead.
"""
struct YeeMesh1D{T<:AbstractFloat}
    ey::Vector{T}
    ez::Vector{T}
    hy::Vector{T}
    hz::Vector{T}
    N::Int
    function YeeMesh1D{T}(N::Integer) where {T}
        new{T}(zeros(T, N+1), zeros(T, N+1), zeros(T, N), zeros(T, N), N)
    end
end

"""
    PML(N, σ_max, Δx, Δt)

Two absorbing layers of `N` cells each, one at either end of the grid, with the
conductivity rising as the cube of the depth to `σ_max` at the wall. `N = 0`
gives no layer at all. `N ≥ 0` and `σ_max ≥ 0` are checked.

The element type is `float` of the arguments' promoted type, so an integer
`σ_max` is accepted: `PML(10, 1000, Δx, Δt)` is `PML(10, 1e3, Δx, Δt)`, to the
bit. Prefer the keyword form, [`PML(; N, σ_max, Δx, Δt)`](@ref PML).
"""
struct PML{T<:AbstractFloat}
    N::Int
    σ_max::T
    r₁::Vector{T}
    r₂::Vector{T}
    function PML(N::Integer, σ_max::Real, Δx::Real, Δt::Real)
        N ≥ 0 || throw(ArgumentError("PML needs N ≥ 0 cells, got $N"))
        σ_max ≥ 0 || throw(ArgumentError("PML needs σ_max ≥ 0, got $σ_max"))
        T = float(promote_type(typeof(σ_max), typeof(Δx), typeof(Δt)))
        # On the arguments as given, as before the element type was computed:
        # `new` converts what comes out, which for Float64 arguments is nothing.
        σ = [σ_max*(i/2N)^3 for i = 1:2N]
        r₁ = exp.(-Δt.*σ)
        # (1 - exp(-Δtσ))/(Δxσ), through `expm1` so that it does not cancel as
        # Δtσ → 0, and at σ = 0 its limit Δt/Δx, the interior coefficient.
        r₂ = [iszero(s) ? Δt/Δx : -expm1(-Δt*s)/(Δx*s) for s in σ]
        new{T}(N, σ_max, r₁, r₂)
    end
    # The same layer with its constants rounded to `T`, for an operator whose
    # step works in `T`. From the coefficients, not from `σ_max`, so that a
    # layer computed in a wider type is rounded once.
    PML{T}(p::PML) where {T<:AbstractFloat} = new{T}(p.N, p.σ_max, p.r₁, p.r₂)
end

"""
    PML(; N, σ_max, Δx, Δt)

Keyword form of the [`PML`](@ref) constructor. The positional form takes
`Δx` before `Δt`, and the two are trivially swappable at a call site --
the solver's own default once had them the wrong way round, from the day it
was written. Prefer this form.
"""
PML(; N, σ_max, Δx, Δt) = PML(N, σ_max, Δx, Δt)

"""
    Yee1D(; Δx, Δt, cfl = Δt/Δx, source, x_min = 0.0,
            pml = PML(; N = 10, σ_max = 1e3, Δx, Δt))

The Yee update on a uniform grid of spacing `Δx`, stepped by `Δt`, for
[`advance!`](@ref). Its fields are those keywords. The step works in one
floating-point type `T`, `float` of the promotion of `Δx` and `Δt`: `cfl`, `Δx`,
`Δt` and the coefficients of `pml` are converted to it, so `Float32` steps give a
`Float32` operator whatever type the literal `cfl` or the default layer was
written in. `x_min` is only a coordinate for the source and keeps its own
floating-point type, so the `x` the source sees is computed as before.

`cfl` must equal `Δt/Δx`: the interior update uses `cfl` and the absorbing
layer `Δt/Δx`, so two different values would put an impedance step at the edge
of the layer. It is checked in `T`, to a relative `max(1e-12, 4eps(T))`, by
every constructor, the positional `Yee1D(cfl, Δx, Δt, x_min, pml, source)`
included. It is a keyword of its own because `cfl*Δx/Δx` need not round back
to `cfl`, and a caller that sets the Courant number gets exactly that number,
rounded to `T`, in the interior.

`source` is a named tuple `(y, z)` of functions `(t, x) -> amplitude`,
injected one-way (rightwards) at the first interior node `pml.N + 2`. Their
`x` argument is `x_min + Δx` for the electric field and `x_min + 1.5Δx` for the
magnetic one, so `x_min` is the coordinate the source is measured from. A
source that is not wanted is `(t, x) -> 0.0`.

The value holds no fields and no scratch, so it can be shared by any number of
meshes and tasks; [`workspace`](@ref Vasilek.workspace) returns `nothing` for
it.
"""
struct Yee1D{T<:AbstractFloat, X<:AbstractFloat, S}
    cfl::T
    Δx::T
    Δt::T
    x_min::X
    pml::PML{T}
    source::S
    function Yee1D(cfl::Real, Δx::Real, Δt::Real, x_min::Real, pml::PML, source)
        T = float(promote_type(typeof(Δx), typeof(Δt)))
        c, h, τ = T(cfl), T(Δx), T(Δt)
        isapprox(c, τ/h; rtol = max(1e-12, 4eps(T))) || throw(ArgumentError(
            "cfl = $cfl but Δt/Δx = $(τ/h): the interior and the absorbing " *
            "layer would use different Courant numbers"))
        x₀ = float(x_min)
        return new{T, typeof(x₀), typeof(source)}(c, h, τ, x₀, PML{T}(pml), source)
    end
end

Yee1D(; Δx, Δt, cfl = Δt/Δx, source, x_min = 0.0,
        pml = PML(; N = 10, σ_max = 1e3, Δx = Δx, Δt = Δt)) =
    Yee1D(cfl, Δx, Δt, x_min, pml, source)

"""
    workspace(op::Yee1D, args...)

`nothing`: the Yee update needs no scratch. Any further arguments are accepted
and ignored, so generic code can call `workspace(op, n)`.
"""
workspace(::Yee1D, args...) = nothing

@noinline _err_fit(Nx, NP) = throw(ArgumentError(
    "a mesh of $Nx cells cannot hold two absorbing layers of $NP cells " *
    "and an interior; need N ≥ $(2*NP + 2)"))

@noinline _err_current(Nx, ny, nz) = throw(DimensionMismatch(
    "a mesh of $Nx cells needs currents on its $(Nx + 1) electric nodes; " *
    "got length(j.y) = $ny and length(j.z) = $nz"))

"""
    advance!(mesh::YeeMesh1D, op::Yee1D, t, j)

One leapfrog step of `mesh` under `op`, ending at time `t`, and `mesh` back:
the source of `op` is injected at `t` (and at `t + Δt/2` for the magnetic
field), then `E` is advanced with the current and `H` after it.

`j` is a named tuple `(y, z)` of arrays indexed like `mesh.ey`, and is added
**straight into the field**, so the caller owes it the time step: the argument
is `-J Δt`, not `J`. See `docs/src/normalization.md`.

Only the interior nodes `2:N` are driven. The two end nodes are PEC boundaries
(see [`YeeMesh1D`](@ref)), so `j.y[1]`, `j.z[1]` and the last entry of each are
ignored -- a current cannot be injected into a perfect conductor. The arrays
still span the full `N+1` nodes so that they index alongside `mesh.ey`.

The mesh must hold both absorbing layers of `op` and an interior,
`mesh.N ≥ 2*op.pml.N + 2`, and `j.y` and `j.z` must have the `N+1` entries
of `mesh.ey`; both are checked here, on every call and before anything is
written, since the operator is built without a mesh.
"""
function advance!(mesh::YeeMesh1D, op::Yee1D, t, j)
    mesh.N ≥ 2*op.pml.N + 2 || _err_fit(mesh.N, op.pml.N)
    length(j.y) == length(j.z) == mesh.N + 1 || _err_current(mesh.N, length(j.y), length(j.z))
    _inject!(mesh, op, t)
    cfl = op.cfl; pml = op.pml; Nx = mesh.N
    _update_e!(mesh.ey, mesh.hz, j.y, -, cfl, pml, Nx)
    _update_e!(mesh.ez, mesh.hy, j.z, +, cfl, pml, Nx)
    _update_h!(mesh.hz, mesh.ey, -, cfl, pml, Nx)
    _update_h!(mesh.hy, mesh.ez, +, cfl, pml, Nx)
    return mesh
end

# The steps of `advance!`. The source is injected first, then each
# polarisation's electric field and, after both, each one's magnetic field.

function _inject!(mesh, op, t)
    f = mesh; cfl = op.cfl; Δt = op.Δt; Δx = op.Δx; x_min = op.x_min
    pml = op.pml; pulse_shape = op.source
    f.ey[pml.N+2] -= cfl*pulse_shape.y(t, x_min + Δx)
    f.ez[pml.N+2] -= cfl*pulse_shape.z(t, x_min + Δx)

    # A right-going wave has hz = ey but hy = -ez: the (ez, hy) pair is
    # (ey, -hz). With `-=` here the z source launched its pulse to the left,
    # into the absorbing layer.
    f.hz[pml.N+2] -= cfl*pulse_shape.y(t + 0.5*Δt, x_min + 1.5*Δx)
    f.hy[pml.N+2] += cfl*pulse_shape.z(t + 0.5*Δt, x_min + 1.5*Δx)
    return nothing
end

# One polarisation's update: `e` on the nodes, `h` on the cells, and `±` the
# sign the curl enters with -- `-` for (ey, hz) and `+` for (ez, hy), the pair
# (ez, hy) being (ey, -hz). `±` is `+` or `-` itself, not a factor of ±1, so
# `e ± cfl*d` is the same operation the two hand-written copies performed.

# `e` on the interior nodes 2:Nx, then the current. The end nodes are PEC.
function _update_e!(e, h, j, ±, cfl, pml, Nx)
    for i = 2:pml.N+1
        k = 1 + 2*(pml.N-i+1)
        e[i] = pml.r₁[k]*e[i] ± pml.r₂[k]*(h[i] - h[i-1])
    end
    for i = pml.N+2:Nx-pml.N
        e[i] = e[i] ± cfl*(h[i] - h[i-1])
    end
    for i = Nx-pml.N+1:Nx
        k = 1 + 2*(i-Nx+pml.N-1)
        e[i] = pml.r₁[k]*e[i] ± pml.r₂[k]*(h[i] - h[i-1])
    end
    for i = 2:Nx
        e[i] += j[i]
    end
    return nothing
end

# `h` on all Nx cells. The last reads `e[end]`, the right-hand PEC node.
function _update_h!(h, e, ±, cfl, pml, Nx)
    for i = 1:pml.N
        k = 2*(pml.N-i+1)
        h[i] = pml.r₁[k]*h[i] ± pml.r₂[k]*(e[i+1] - e[i])
    end
    for i = pml.N+1:Nx-pml.N
        h[i] = h[i] ± cfl*(e[i+1] - e[i])
    end
    for i = Nx-pml.N+1:Nx
        k = 2*(i-Nx+pml.N)
        h[i] = pml.r₁[k]*h[i] ± pml.r₂[k]*(e[i+1] - e[i])
    end
    return nothing
end

end
