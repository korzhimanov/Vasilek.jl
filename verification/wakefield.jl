# Laser wakefield excitation in a 1D1V plasma slab.
#
#     julia --project=verification verification/wakefield.jl
#
# Writes wakefield-*.png beside this script.
#
# The physics lives in `wakefield` in `test/verification_harness.jl`, so that
# this script and the test asserting its claims run the same code rather than
# two copies that drift. Before that extraction there was nothing for a test to
# call, which is why the README's "runs and is stable" went unasserted for as
# long as it did.
#
# KNOWN INCOMPLETE: there is no ponderomotive coupling. The laser never enters
# the longitudinal push, so the wake below is the slab edges relaxing rather
# than a laser-driven wave, and `peak wake field` and `Δε/ε` come out
# bit-identical whether the transverse current is right, wrong by Δt, or wrong
# by thirty-two orders of magnitude. Closing that needs the ponderomotive force
# −∇(pʸ² + pᶻ²)/2γ in the momentum advection, which is a modelling decision for
# the author. See the docstring on `wakefield` for the rest.

using Plots

include(joinpath(@__DIR__, "..", "test", "verification_harness.jl"))

r = wakefield()

here = @__DIR__

# The laser field over space and time. This used to be `heatmap(x, t, em.ey)`,
# passing the final snapshot -- a length-Nx vector -- where a Nt×Nx matrix was
# wanted. It never raised, because a script never renders the plot.
savefig(heatmap(r.x/2π, r.t/2π, r.ey; xlabel = "x/2π", ylabel = "t/2π",
                title = "eʸ (laser)"),
        joinpath(here, "wakefield-laser.png"))

savefig(heatmap(r.x/2π, r.t/2π, r.ex; xlabel = "x/2π", ylabel = "t/2π",
                title = "eˣ (wake)"),
        joinpath(here, "wakefield-wake.png"))

savefig(heatmap(r.x/2π, r.t/2π, r.n; xlabel = "x/2π", ylabel = "t/2π",
                title = "nₑ"),
        joinpath(here, "wakefield-density.png"))

# Energy diagnostics were accumulated but never shown.
energy = plot(r.t/2π, (r.ε .- r.ε[1])./r.ε[1]; label = "Δε/ε", xlabel = "t/2π")
plot!(energy, r.t/2π, r.ε_e./r.ε[1]; label = "electrostatic / ε₀")
savefig(energy, joinpath(here, "wakefield-energy.png"))

println("final Δε/ε      = ", (r.ε[end] - r.ε[1])/r.ε[1])
println("peak wake field = ", maximum(abs, r.ex))
println("peak laser field= ", maximum(abs, r.ey))
println("wrote wakefield-{laser,wake,density,energy}.png to ", here)
