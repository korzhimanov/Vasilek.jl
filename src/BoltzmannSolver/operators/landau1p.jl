"""
    ∂f∂v(f, v, k)

Derivative of `f` with respect to `v` at index `k`: centred in the interior,
one-sided at either end.

Factored out deliberately. This used to be written twice inside the solver, and
the second copy branched on the outer loop index while claiming to
differentiate at the inner one, so the two "different" derivatives were always
equal and the collision operator computed something else entirely.
"""
function ∂f∂v(f, v, k)
    if k == 1
        return (f[2] - f[1])/(v[2] - v[1])
    elseif k == length(f)
        return (f[end] - f[end-1])/(v[end] - v[end-1])
    else
        return (f[k+1] - f[k-1])/(v[k+1] - v[k-1])
    end
end

"""
    Landau1P(A; L = 20.0, Tₜ = 1e-3)

Landau collision integral in one-dimensional velocity space.
`A = 4πe⁴N₀/(m²v₀³ω)` is a normalizing constant and `L` the Coulomb logarithm.
`Tₜ` is a transversal temperature used to estimate `|vₛ-vₚ|² - (vₛₓ-vₚₓ)²` as
`≈2Tₜ`.

!!! warning "Experimental"
    The closure is inconsistent: the numerator carries the transversal estimate
    `2Tₜ` while the denominator uses the purely longitudinal `|vᵢ-vⱼ|³`, leaving
    a non-integrable singularity at `i ≈ j`. The collision integral grows under
    refinement of the velocity grid rather than converging, and a Maxwellian is
    not a stationary point. A consistent closure would use `(Δv² + 2Tₜ)^(3/2)`.

    The update also differences a cell-centred `I` rather than staggered fluxes,
    so mass is not conserved to machine precision even in principle.

    Not exported.
"""
struct Landau1P{T<:AbstractFloat} <: AbstractCollisionOperator
    A::T
    L::T
    Tₜ::T
end
Landau1P(A; L = 20.0, Tₜ = 1e-3) = Landau1P(promote(float(A), float(L), float(Tₜ))...)

struct Landau1PWorkspace{T}
    I::Vector{T}
    J::Vector{T}
end
workspace(::Landau1P, n::Integer) = Landau1PWorkspace(Vector{Float64}(undef, n), Vector{Float64}(undef, n))

function collide!(dest, src, op::Landau1P, v, Δt, ws::Landau1PWorkspace)
    I, J = ws.I, ws.J
    n = length(src)

    for i in eachindex(src)
        Δfᵢ = ∂f∂v(src, v, i)
        for j in eachindex(src)
            if i == j
                J[j] = zero(eltype(J))   # the kernel is singular here
            else
                J[j] = (src[i]*∂f∂v(src, v, j) - src[j]*Δfᵢ)*2*op.Tₜ/abs(v[i]-v[j])^3
            end
        end
        I[i] = op.L*op.A*integrate(v, J)
    end

    for i in eachindex(src)
        if i == 1
            dest[1] = src[1] - Δt*(I[2] - I[1])/(v[2] - v[1])
        elseif i == n
            dest[end] = src[end] - Δt*(I[end] - I[end-1])/(v[end] - v[end-1])
        else
            dest[i] = src[i] - Δt*(I[i+1] - I[i-1])/(v[i+1] - v[i-1])
        end
    end
    return dest
end
