# Laser wakefield excitation in a 1D1V plasma slab.
#
#     julia --project=verification verification/wakefield.jl
#
# Writes wakefield-*.png beside this script.

using Plots
using NumericalIntegration

using Vasilek
using Vasilek: FDTD1D, PoissonFourier1D

Δx = 0.1*2π
Δt = 0.05*Δx
Δp = 0.1

x_min = -5.0*2π
box_length = 20.0*2π
plasma_thickness = 10.0*2π
plasma_temperature = 0.1
plasma_density = 0.1

laser_amplitude = 1.0
laser_duration = 5*2π

total_time = 2π*10

x = collect(x_min:Δx:(box_length + x_min))
p = collect(-laser_amplitude*4:Δp:laser_amplitude*4)

Nx = length(x)
Np = length(p)

f = plasma_density/sqrt(2π*plasma_temperature) *
    (@. exp(-0.5*(p)^2/plasma_temperature)) *
    (@. 0.5*(tanh(x) - tanh(x - plasma_thickness)))'

# Immobile neutralising ions, at the initial electron density.
nᵢ = integrate(p, f)

g = similar(f)

# PFC holds no arrays, so one value serves every line of both sweeps.
advection = PFC(fmin = 0.0, fmax = maximum(f))

em = FDTD1D.YeeMesh1D{Float64}(Nx - 1)

pulse_shape = (
    y = (t, x) -> exp(-((x - t)/laser_duration)^2)*sin(x - t),
    z = (t, x) -> 0.0
)

advance_fields! = FDTD1D.make_advance_fields(em, Δt/Δx, pulse_shape, Δt, Δx, x_min,
                                             FDTD1D.PML(; N = 0, σ_max = 1.0, Δx = Δx, Δt = Δt))

t = collect(0.0:Δt:total_time)
Nt = length(t)

n = zeros(Nt, Nx)          # electron density history
ey = zeros(Nt, Nx)         # transverse laser field history
ex = zeros(Nt, Nx)         # longitudinal wake history
ε_e = zeros(Nt)            # electrostatic energy
ε = zeros(Nt)              # total energy

solve_poisson! = PoissonFourier1D.generate_solver(nᵢ, Δx)

e = similar(x)
n[1, :] = nᵢ
solve_poisson!(e, n[1, :] - nᵢ)
ex[1, :] = e
ey[1, :] = em.ey
ε_e[1] = integrate(x, e.^2)
ε[1] = integrate(x, integrate(p, @. f*p^2)) + ε_e[1]

pʸ = zeros(Nx)
pᶻ = zeros(Nx)

for k in 2:Nt
    for j = 1:Np
        advect!(view(g, j, :), view(f, j, :), advection, p[j]*Δt/Δx)
    end

    n[k, :] = integrate(p, g)
    solve_poisson!(e, n[k, :] - nᵢ)

    pʸ .= pʸ .+ em.ey.*Δt
    pᶻ .= pᶻ .+ em.ez.*Δt

    # Transverse current. Two corrections against what this line used to be,
    # `(y = pʸ.*n[k, :], z = pᶻ.*n[k, :])`:
    #
    #   * Δt, because make_advance_fields adds the argument straight into ey
    #     (`f.ey[i] += jy[i]`), so the caller owes it the time step. Without it
    #     the peak field reached 1.0e22.
    #   * the sign, so that ∂p/∂t = +e and ∂e/∂t = -n·p together give an
    #     oscillation rather than exponential growth. That pairing is the one
    #     the longitudinal direction already uses, and the plasma-oscillation
    #     notebook verifies it reproduces the right plasma frequency -- so it
    #     is taken from working code rather than chosen here. With the sign as
    #     written the peak field was still 44.
    #
    # With both, the peak transverse field is 0.38 and bounded.
    #
    # KNOWN INCOMPLETE: there is no ponderomotive coupling. The laser never
    # enters the longitudinal push, so the wake below is the slab edges
    # relaxing, not a laser-driven wave -- `peak wake field` and `Δε/ε` come
    # out bit-identical whether the transverse current is right, wrong by Δt,
    # or wrong by 32 orders of magnitude. Closing that needs the ponderomotive
    # force -∇(pʸ² + pᶻ²)/2γ in the momentum advection, which is a modelling
    # decision for the author.
    #
    # NOTE also that the current is taken through momentum rather than
    # velocity. At laser_amplitude = 1.0 the motion is relativistic (a₀ = 1)
    # and 1/γ is not close to unity, with γ = sqrt(1 + pₓ² + pʸ² + pᶻ²).
    advance_fields!(k*Δt, (y = -pʸ.*n[k, :].*Δt, z = -pᶻ.*n[k, :].*Δt))

    for i = 1:Nx
        advect!(view(f, :, i), view(g, :, i), advection, e[i]*Δt/Δp)
    end

    ey[k, :] = em.ey
    ex[k, :] = e
    ε_e[k] = integrate(x, e.^2)
    ε[k] = integrate(x, integrate(p, @. f*p^2)) + ε_e[k]
end

here = @__DIR__

# The laser field over space and time. This used to be `heatmap(x, t, em.ey)`,
# passing the final snapshot -- a length-Nx vector -- where a Nt×Nx matrix was
# wanted. It never raised, because a script never renders the plot.
savefig(heatmap(x/2π, t/2π, ey; xlabel = "x/2π", ylabel = "t/2π", title = "eʸ (laser)"),
        joinpath(here, "wakefield-laser.png"))

savefig(heatmap(x/2π, t/2π, ex; xlabel = "x/2π", ylabel = "t/2π", title = "eˣ (wake)"),
        joinpath(here, "wakefield-wake.png"))

savefig(heatmap(x/2π, t/2π, n; xlabel = "x/2π", ylabel = "t/2π", title = "nₑ"),
        joinpath(here, "wakefield-density.png"))

# Energy diagnostics were accumulated but never shown.
energy = plot(t/2π, (ε .- ε[1])./ε[1]; label = "Δε/ε", xlabel = "t/2π")
plot!(energy, t/2π, ε_e./ε[1]; label = "electrostatic / ε₀")
savefig(energy, joinpath(here, "wakefield-energy.png"))

println("final Δε/ε      = ", (ε[end] - ε[1])/ε[1])
println("peak wake field = ", maximum(abs, ex))
println("peak laser field= ", maximum(abs, ey))
println("wrote wakefield-{laser,wake,density,energy}.png to ", here)
