module PoissonFourier1D
export solve_poisson!

using FFTW

function solve_poisson!(e, ω, ρ, Δx)
    F = FFTW.rfft(ρ)
    φ = FFTW.irfft(F./(-ω.^2), length(ρ))
    e[:] = vcat(0.5*[φ[2]-φ[end]], 0.5*(φ[3:end] - φ[1:end-2]), 0.5*[φ[1]-φ[end-1]])./Δx
end

end
