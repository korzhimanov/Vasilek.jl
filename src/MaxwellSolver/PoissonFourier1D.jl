module PoissonFourier1D

using FFTW, LinearAlgebra

function generate_solver(ρ₀, Δx)
    PF = FFTW.plan_rfft(ρ₀)
    F = FFTW.rfft(ρ₀)
    ω = collect(2π*rfftfreq(length(ρ₀)))
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
        for i = 2:length(φ)-1
            e[i] = 0.5*(φ[i+1] - φ[i-1])/Δx
        end
        e[end] = 0.5*(φ[1]-φ[end-1])/Δx
        return e
    end
    return solve!
end

end
