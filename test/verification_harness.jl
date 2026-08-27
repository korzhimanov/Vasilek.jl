# Shared 1D1V Vlasov–Poisson driver for the extended verification tests.
# Deliberately close to what the verification notebooks do, so the two agree.

using Vasilek
using Vasilek: StrangSplitting
using NumericalIntegration, FFTW

"Spectral Poisson solve on a uniform x grid, e = -dφ/dx with φ'' = -ρ."
function make_poisson(x)
    Δx = x[2] - x[1]
    L = x[end] - x[1] + Δx
    ω = 2π*collect(0.0:1.0/L:0.5/Δx)
    ω[1] = ω[2]
    return function (e, ρ)
        F = FFTW.rfft(ρ)
        φ = FFTW.irfft(F ./ (-ω.^2), length(ρ))
        e[:] = vcat(0.5*[φ[2] - φ[end]],
                    0.5*(φ[3:end] - φ[1:end-2]),
                    0.5*[φ[1] - φ[end-1]]) ./ Δx
        return e
    end
end

cell_widths(z) = vcat([z[2]-z[1]], 0.5*(z[3:end] - z[1:end-2]), [z[end]-z[end-1]])

"""
    vlasov_poisson(x, v, f₀, t)

Strang-split Vlasov–Poisson on a static grid, returning the electric energy
history and the total energy history.
"""
function vlasov_poisson(x, v, f₀, t)
    Δx = cell_widths(x)
    Δv = cell_widths(v)
    scheme_x = PFCNonUniform(Δx; fmin = 0.0, fmax = 1.0)
    scheme_v = PFCNonUniform(Δv; fmin = 0.0, fmax = 1.0)
    ws_x = workspace(scheme_x, length(Δx))
    ws_v = workspace(scheme_v, length(Δv))
    buf_x = Vector{Float64}(undef, length(Δx))
    buf_v = Vector{Float64}(undef, length(Δv))

    advect_x! = (col, α) -> (advect!(buf_x, col, scheme_x, α, ws_x); copyto!(col, buf_x))
    advect_v! = (col, α) -> (advect!(buf_v, col, scheme_v, α, ws_v); copyto!(col, buf_v))

    fᵢ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. Δx/Δx)'
    nᵢ = integrate(v, fᵢ)
    Nᵢ = integrate(x, nᵢ)

    f = copy(f₀)
    f .*= Nᵢ/integrate(x, integrate(v, f))
    g = f'

    solve_poisson! = make_poisson(x)
    e = similar(x)
    ε_e = similar(t)
    ε = similar(t)

    for k in 1:length(t)-1
        Δt = t[k+1] - t[k]
        vΔt(_) = v*Δt
        function eΔt(ff)
            nₖ = vec(sum(ff'.*Δv, dims = 1))
            solve_poisson!(e, nₖ - nᵢ)
            return e*Δt
        end
        StrangSplitting.make_time_step_2d!((g, f), (vΔt, eΔt), (advect_x!, advect_v!))
        ε_e[k] = integrate(x, e.^2)
        ε[k] = integrate(x, integrate(v, @. f*v^2)) + ε_e[k]
    end
    ε_e[end] = ε_e[end-1]
    ε[end] = ε[end-1]
    return ε_e, ε
end
