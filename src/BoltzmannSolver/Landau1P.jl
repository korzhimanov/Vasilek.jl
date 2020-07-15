module Landau1P

using NumericalIntegration

"""
    generate_solver(f₀, f, v_, Δt, A)

Generates solver for Landau integral in 1-dimensional velocity space. ``A = 4πe⁴N₀/(m²v₀³ω)`` is a normalizing constant.
"""
function generate_solver(f₀, f, v_, Δt, A)
    v = copy(v_)
    L = 20.0 # Coulomb logarithm
    Tₜ = 1e-3 # transversal temperature needed to estimate |vₛ-vₚ|² - (vₛₓ-vₚₓ)² as ≈2Tₜ

    I = similar(f)
    J = similar(v)

    function solve!()
        for i in eachindex(f)
            if i == 1
                Δfᵢ = (f₀[2]-f₀[1])/(v[2]-v[1])
            elseif i == length(f)
                Δfᵢ = (f₀[end]-f₀[end-1])/(v[end]-v[end-1])
            else
                Δfᵢ = (f₀[i+1]-f₀[i-1])/(v[i+1]-v[i-1])
            end
            for j in eachindex(f₀)
                if i != j
                    if i == 1
                        Δfⱼ = (f₀[2]-f₀[1])/(v[2]-v[1])
                    elseif i == length(f)
                        Δfⱼ = (f₀[end]-f₀[end-1])/(v[end]-v[end-1])
                    else
                        Δfⱼ = (f₀[i+1]-f₀[i-1])/(v[i+1]-v[i-1])
                    end
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
