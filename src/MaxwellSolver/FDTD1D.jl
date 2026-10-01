module FDTD1D
export make_advance_fields, YeeMesh1D

"""
    YeeMesh1D{T}(N)

Staggered fields on `N` cells: `N+1` electric nodes with `N` magnetic cells
between them, all zeroed.

The two end nodes, `ey[1]`/`ez[1]` and `ey[end]`/`ez[end]`, are **perfect
electric conductor boundaries**. No update touches them — the interior loop
runs `pml.N+2 : Nx-pml.N` and the two absorbing-layer loops stop short of both
ends — so they hold the zero they are constructed with, for all time. They are
not dead storage: `update_hz!` reads `ey[end]`, which is how the boundary
condition enters the solution.

Writing a nonzero value into either end node therefore does not seed a wave,
it silently changes the boundary condition to `E = const`. Seed the interior
instead.
"""
struct YeeMesh1D{T<:AbstractFloat, S<:Integer}
    ey::Vector{T}
    ez::Vector{T}
    hy::Vector{T}
    hz::Vector{T}
    N::S
    function YeeMesh1D{T}(N::S) where {T, S}
        new{T,S}(zeros(T, N+1), zeros(T, N+1), zeros(T, N), zeros(T, N), N)
    end
end

struct PML{T<:Integer, S<:AbstractFloat}
    N::T
    σ_max::S
    r₁::Vector{S}
    r₂::Vector{S}
    function PML(N::T, σ_max::S, Δx, Δt) where {T, S}
        N ≥ 0 || throw(ArgumentError("PML needs N ≥ 0 cells, got $N"))
        σ_max ≥ 0 || throw(ArgumentError("PML needs σ_max ≥ 0, got $σ_max"))
        σ = [σ_max*(i/2N)^3 for i = 1:2N]
        r₁ = exp.(-Δt.*σ)
        # (1 - exp(-Δtσ))/(Δxσ), through `expm1` so that it does not cancel as
        # Δtσ → 0, and at σ = 0 its limit Δt/Δx, the interior coefficient.
        r₂ = [iszero(s) ? Δt/Δx : -expm1(-Δt*s)/(Δx*s) for s in σ]
        new{T,S}(N, σ_max, r₁, r₂)
    end
end

"""
    PML(; N, σ_max, Δx, Δt)

Keyword form of the [`PML`](@ref) constructor. The positional form takes
`Δx` before `Δt`, and the two are trivially swappable at a call site --
the default argument of `make_advance_fields` had them the wrong way
round from the day it was written. Prefer this form.
"""
PML(; N, σ_max, Δx, Δt) = PML(N, σ_max, Δx, Δt)

"""
    make_advance_fields(f, cfl, pulse_shape, Δt, Δx, x_min[, pml])

Build `advance_fields!(t, j)`, one leapfrog step on the mesh `f`.

`cfl` must equal `Δt/Δx`: the interior update uses `cfl` and the absorbing
layer `Δt/Δx`, so two different values would put an impedance step at the edge
of the layer. It is checked.

`pulse_shape` is a named tuple `(y, z)` of functions `(t, x) -> amplitude`,
injected one-way (rightwards) at the first interior node `pml.N + 2`. Their
`x` argument is `x_min + Δx` for the electric field and `x_min + 1.5Δx` for the
magnetic one, so `x_min` is the coordinate the source is measured from.

`j` is a named tuple `(y, z)` of arrays indexed like `f.ey`, and is added
**straight into the field**, so the caller owes it the time step: the argument
is `-J Δt`, not `J`. See `docs/normalization.md`.

Only the interior nodes `2:N` are driven. The two end nodes are PEC boundaries
(see [`YeeMesh1D`](@ref)), so `j.y[1]`, `j.z[1]` and the last entry of each are
ignored -- a current cannot be injected into a perfect conductor. The arrays
still span the full `N+1` nodes so that they index alongside `f.ey`.
"""
function make_advance_fields(f::YeeMesh1D{T,S}, cfl, pulse_shape, Δt, Δx, x_min, pml::PML = PML(; N = 10, σ_max = 1e3, Δx = Δx, Δt = Δt)) where {T,S}
    Nx = f.N
    isapprox(cfl, Δt/Δx; rtol = 1e-12) || throw(ArgumentError(
        "cfl = $cfl but Δt/Δx = $(Δt/Δx): the interior and the absorbing layer " *
        "would use different Courant numbers"))
    Nx ≥ 2*pml.N + 2 || throw(ArgumentError(
        "a mesh of $Nx cells cannot hold two absorbing layers of $(pml.N) cells " *
        "and an interior; need N ≥ $(2*pml.N + 2)"))

    function generate_fields_x_min!(t)
        f.ey[pml.N+2] -= cfl*pulse_shape.y(t, x_min + Δx)
        f.ez[pml.N+2] -= cfl*pulse_shape.z(t, x_min + Δx)
        
        # A right-going wave has hz = ey but hy = -ez: the (ez, hy) pair is
        # (ey, -hz). With `-=` here the z source launched its pulse to the left,
        # into the absorbing layer.
        f.hz[pml.N+2] -= cfl*pulse_shape.y(t + 0.5*Δt, x_min + 1.5*Δx)
        f.hy[pml.N+2] += cfl*pulse_shape.z(t + 0.5*Δt, x_min + 1.5*Δx)
    end

    function update_ey!(jy)
        for i = 2:pml.N+1
            f.ey[i] = pml.r₁[1+2*(pml.N-i+1)]*f.ey[i] - pml.r₂[1+2*(pml.N-i+1)]*(f.hz[i] - f.hz[i-1])
        end
        for i = pml.N+2:Nx-pml.N
            f.ey[i] -= cfl*(f.hz[i] - f.hz[i-1])
        end
        for i = Nx-pml.N+1:Nx
            f.ey[i] = pml.r₁[1+2*(i-Nx+pml.N-1)]*f.ey[i] - pml.r₂[1+2*(i-Nx+pml.N-1)]*(f.hz[i] - f.hz[i-1])
        end
        for i = 2:Nx
            f.ey[i] += jy[i]
        end
    end
    
    function update_ez!(jz)
        for i = 2:pml.N+1
            f.ez[i] = pml.r₁[1+2*(pml.N-i+1)]*f.ez[i] + pml.r₂[1+2*(pml.N-i+1)]*(f.hy[i] - f.hy[i-1])
        end
        for i = pml.N+2:Nx-pml.N
            f.ez[i] += cfl*(f.hy[i] - f.hy[i-1])
        end
        for i = Nx-pml.N+1:Nx
            f.ez[i] = pml.r₁[1+2*(i-Nx+pml.N-1)]*f.ez[i] + pml.r₂[1+2*(i-Nx+pml.N-1)]*(f.hy[i] - f.hy[i-1])
        end
        for i = 2:Nx
            f.ez[i] += jz[i]
        end
    end
    
    function update_hy!()
        for i = 1:pml.N
            f.hy[i] = pml.r₁[2*(pml.N-i+1)]*f.hy[i] + pml.r₂[2*(pml.N-i+1)]*(f.ez[i+1] - f.ez[i])
        end
        for i = pml.N+1:Nx-pml.N
            f.hy[i] += cfl*(f.ez[i+1] - f.ez[i])
        end
        for i = Nx-pml.N+1:Nx
            f.hy[i] = pml.r₁[2*(i-Nx+pml.N)]*f.hy[i] + pml.r₂[2*(i-Nx+pml.N)]*(f.ez[i+1] - f.ez[i])
        end
    end
    
    function update_hz!()
        for i = 1:pml.N
            f.hz[i] = pml.r₁[2*(pml.N-i+1)]*f.hz[i] - pml.r₂[2*(pml.N-i+1)]*(f.ey[i+1] - f.ey[i])
        end
        for i = pml.N+1:Nx-pml.N
            f.hz[i] -= cfl*(f.ey[i+1] - f.ey[i])
        end
        for i = Nx-pml.N+1:Nx
            f.hz[i] = pml.r₁[2*(i-Nx+pml.N)]*f.hz[i] - pml.r₂[2*(i-Nx+pml.N)]*(f.ey[i+1] - f.ey[i])
        end
    end

    function make_step!(j)
        update_ey!(j.y)
        update_ez!(j.z)
        
        update_hz!()
        update_hy!()
    end
    
    function advance_fields!(t, j)
        generate_fields_x_min!(t)
        make_step!(j)
    end
    
    return advance_fields!
end

end
