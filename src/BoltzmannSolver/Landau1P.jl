module Landau1P

using NumericalIntegration

"""
    ∂f∂v(f, v, k)

Derivative of `f` with respect to `v` at index `k`: centred in the interior,
one-sided at either end.

Factored out deliberately. This expression used to be written twice inside
`solve!`, and the second copy branched on the outer loop index while claiming
to differentiate at the inner one -- so the two "different" derivatives were
always equal and the collision operator computed something else entirely.
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
    generate_solver(f₀, f, v_, Δt, A; L = 20.0, Tₜ = 1e-3)

Generates a solver for the Landau integral in 1-dimensional velocity space.
``A = 4πe⁴N₀/(m²v₀³ω)`` is a normalizing constant, `L` the Coulomb logarithm.

`Tₜ` is a transversal temperature used to estimate ``|vₛ-vₚ|² - (vₛₓ-vₚₚ)²``
as ``≈2Tₜ``.

!!! warning "Experimental"
    The closure is inconsistent: the numerator carries the transversal
    estimate `2Tₜ` while the denominator uses the purely longitudinal
    `|vᵢ-vⱼ|³`. That leaves a non-integrable singularity at `i ≈ j`, so the
    collision integral grows rather than converges under refinement of the
    velocity grid — see the refinement test, which is marked broken. A
    consistent closure would use `(Δv² + 2Tₜ)^(3/2)` in the denominator.

    The update also differences a cell-centred `I` rather than staggered
    fluxes, so mass is not conserved to machine precision even in principle.
"""
function generate_solver(f₀, f, v_, Δt, A; L = 20.0, Tₜ = 1e-3)
    v = copy(v_)

    I = similar(f)
    J = similar(v)

    function solve!()
        for i in eachindex(f)
            Δfᵢ = ∂f∂v(f₀, v, i)
            for j in eachindex(f₀)
                if i == j
                    J[j] = zero(eltype(J))   # the kernel is singular here
                else
                    Δfⱼ = ∂f∂v(f₀, v, j)
                    J[j] = (f₀[i]*Δfⱼ - f₀[j]*Δfᵢ)*2Tₜ/abs(v[i]-v[j])^3
                end
            end
            I[i] = L*A*integrate(v, J)
        end

        for i in eachindex(f)
            if i == 1
                f[1] = f₀[1] - Δt*(I[2] - I[1])/(v[2] - v[1])
            elseif i == length(f)
                f[end] = f₀[end] - Δt*(I[end] - I[end-1])/(v[end] - v[end-1])
            else
                f[i] = f₀[i] - Δt*(I[i+1] - I[i-1])/(v[i+1] - v[i-1])
            end
        end
    end

    return solve!
end
end  # module Landau1P
