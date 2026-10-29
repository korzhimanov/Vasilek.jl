# # A relativistic laser and a thin foil
#
# `VlasovMaxwell1D1P` on the problem it was written for: a plane pulse of
# `a₀ = 20`, 11 fs FWHM at 0.8 µm, on a foil of 40 n_c a few tens of nanometres
# thick, electrons and C⁶⁺ both kinetic, against what can be computed without
# it.
#
#     julia --project=verification verification/laser-foil.jl
#
# Writes laser-foil-*.png beside this script. #src
#
# The setup is `laser_foil` in `test/vm_harness.jl`, which the tests and this
# script share, so that the numbers below and the assertions in
# `test/test_vlasov_maxwell.jl` are about one setup. It is the common 1D setup
# of the transparency study it was built for: the pulse enters at x = −20 µm,
# the foil occupies [0, l], the field is kept on [−20, 15] µm between absorbing
# layers and `f` on [−5, 10] µm, and the transmittance and reflectance are the
# energies of the right- and left-going waves `F = (E_y + B_z)/2`,
# `G = (E_y − B_z)/2` at t = 110.4 fs over the incident energy, measured in a
# vacuum run at 60 fs.

using Plots
using Printf
using Vasilek

include(joinpath(pkgdir(Vasilek), "test", "vm_harness.jl"))

here = @__DIR__;

# ## The pulse in vacuum
#
# The source is `Yee1D`'s one-way injection, which launches minus the function
# it is given; the harness hands it minus the pulse. At 2 nm, 400 cells per
# wavelength and a Courant number of 0.937, the pulse arrives at x = 0 with the
# formula's energy, width and time.

r = vacuum_check()
@printf "W_inc / ∫F²dt of the formula − 1 = %.2e\n" r.W_inc/r.W_formula - 1
@printf "W(t₂, [1, 15] µm) / W_inc − 1     = %.2e\n" r.W_t2/r.W_inc - 1
@printf "leak to the left, ∫G² / W_inc     = %.2e\n" r.G2_t2/r.W_inc
@printf "peak |E_y| at x = 0 / a₀          = %.5f (the formula's sample: %.5f)\n" r.amplitude exp(-(π/2/r.meta.τp)^2)
@printf "FWHM of the intensity at x = 0    = %.4f fs\n" r.fwhm_fs
@printf "peak at x = 0                     = %.3f fs (formula %.3f)\n" r.t_peak_fs 2r.meta.τp*r.meta.u.T + 20/C_UM_PER_FS
@printf "reflection of the right layer     = %.1e\n" r.R_layer

# ## Linear reflection and transmission
#
# At `a₀ = 10⁻³` the foil is a cold dielectric, `ε = 1 − n/n_c`, and its
# reflectance and transmittance are Fresnel's for a layer, averaged over the
# pulse's spectrum (`fresnel_pulse`). With the ions free, a foil of 1 keV
# expands for 85 fs before the pulse's peak arrives, the thinnest by more than
# its own thickness; the reference for those runs is the layer stack of the
# electrons' own profile at that moment (`fresnel_profile`). With ions of
# 10⁹ mₑ the foil keeps its shape to the Debye length, and the reference is the
# uniform layer itself.

rows = []
for l in (10.0, 30.0, 90.0), fixed in (false, true)
    kw = fixed ? (; mi = 1e9, ion_grid = (pmin = -1e-5, pmax = 1e-5, core = 1e-5, Δpmax = 1e-6)) : (;)
    q = fresnel_check(l; kw...)
    push!(rows, (; l, fixed, q.R, q.T, slab = q.slab, profile = q.profile, thin = q.thin))
    @printf "l = %2.0f nm, ions %-5s: R = %.5f T = %.4e | layer R = %.5f T = %.4e | profile R = %.5f T = %.4e | thin foil R = %.4f\n" l (fixed ? "fixed" : "free") q.R q.T q.slab[1] q.slab[2] q.profile[1] q.profile[2] q.thin
end

fig = plot(xlabel = "foil thickness, nm", ylabel = "T", yscale = :log10, legend = :bottomleft)
for fixed in (false, true)
    rs = filter(r -> r.fixed == fixed, rows)
    scatter!(fig, [r.l for r in rs], [r.T for r in rs]; label = fixed ? "run, fixed ions" : "run, free ions")
    plot!(fig, [r.l for r in rs], [(fixed ? r.slab : r.profile)[2] for r in rs];
          label = fixed ? "Fresnel, layer" : "Fresnel, profile at t_f", ls = :dash)
end
savefig(fig, joinpath(here, "laser-foil-fresnel.png"));
#md # ![Linear transmittance against Fresnel](laser-foil-fresnel.png)

# ## The nonlinear run conserves energy
#
# At `a₀ = 20` and ξ = 0.972 the foil turns transparent on the pulse's rising
# edge. The budget between the two probes — the transverse field, the
# longitudinal one, the particles' kinetic energy, what the particles carried
# out of the grid of `f` and what the field carried out of the box — has to
# add up to what came in. At 4 nm, which this script runs, and at the 2 nm of
# the reference run, it does to 0.3% and 0.13%.

n = foil_transmission(; dx_nm = 4.0)
@printf "Δx = 4 nm: T = %.4f R = %.4f A = %.4f, balance %.5f, residual %.2e\n" n.T n.R n.A n.balance n.residual

# The convergence of the transmittance in Δx, Δp, Δt and the scheme takes
# hours, and is the step-10 runner's: `docs/src/verification.md` gives the
# numbers.
println("wrote laser-foil-fresnel.png to ", here) #src
