module PoissonFourier1D
export make_solve_poisson

using FFTW, LinearAlgebra

function generate_solve(ρ₀, Δx)
    PF = FFTW.plan_rfft(ρ₀)
    F = FFTW.rfft(ρ₀)
    ω = collect(LinRange(0.0, π/Δx, length(F)))
    ω[1] = ω[2]
    ω⁻² = -1.0./(ω.^2)
    Pφ = FFTW.plan_irfft(F, length(ρ₀))
    φ = FFTW.irfft(F, length(ρ₀))
    function solve!(e::Vector{T}, ρ::Vector{T}) where {T}
        mul!(F, PF, ρ)
        F[1] = 0.0 + 0.0im
        for i = 2:length(F)
            F[i] = ω⁻²[i] * F[i]
        end
        mul!(φ, Pφ, F)
        e[1] = 0.5*(φ[2]-φ[end])/Δx
        for i = 2:length(e)-1
            e[i] = 0.5*(φ[i+1] - φ[i-1])/Δx
        end
        e[end] = 0.5*(φ[1]-φ[end-1])/Δx
        return e
    end
    return solve!
end

end
