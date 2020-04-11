module PoissonFourier1D
export solve_poisson!

using FFTW

function solve_poisson!(e, ρ, Δx)
    F = FFTW.rfft(ρ)
    ω = collect(LinRange(0.0, π/Δx, length(F)))
    ω[1] = ω[2]
    φ = FFTW.irfft(F./(-ω.^2), length(ρ))
    e[:] = vcat(0.5*[φ[2]-φ[end]], 0.5*(φ[3:end] - φ[1:end-2]), 0.5*[φ[1]-φ[end-1]])./Δx
end

end
