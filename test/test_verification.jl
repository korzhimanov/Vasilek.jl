# Extended verification. About a minute and a half all told -- cheap enough to
# run far more often than nightly, and CI does run it on every PR -- but still
# gated so the default test run stays instant:
#
#     VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
#
# These promote the claims the verification notebooks make in prose into
# assertions. Until they run, those claims are a 2021 HTML file.

include(joinpath(@__DIR__, "verification_harness.jl"))

"""
Roots of the Landau dispersion relation for a Maxwellian, `(γ, ω_r)` by `kλ_D`.

These are the tabulated values, not the `π/(8√2)u³exp(-u²/2)` asymptotic the
notebooks quote -- that gives 0.1145 at k = 0.5 and is not accurate there, and
measuring against it once suggested a 31% error where there is none.
"""
const LANDAU_ROOTS = Dict(0.3 => (0.01262, 1.15985),
                          0.4 => (0.06613, 1.28506),
                          0.5 => (0.15336, 1.41566))

"""
    landau_case(k, Nx, vmax, Δt, tmax)

A single-mode Landau run: two wavelengths of `k` across the box, a Maxwellian
perturbed by 1%, returning `(t, ε_e)`.

**The grid is not free.** The velocity window has to contain the resonance at
`v = ω_r/k` -- 2.83, 3.21 and 3.87 for the three cases below -- because that is
where the damping comes from; and the window then fixes the time step, since the
fastest row runs at `max|v|·Δt/Δx` and `PFCNonUniform` is a finite-volume scheme
with a Courant limit of 1. Widening the window to reach a resonance therefore
costs a smaller `Δt`, not just more velocity points. Every case below runs
between 0.64 and 0.82.
"""
function landau_case(k, Nx, vmax, Δt, tmax)
    L = 2*(2π/k)
    Δx = L/Nx
    x = collect(Δx:Δx:L)
    v = collect(-vmax:0.1:vmax)
    t = collect(0.0:Δt:tmax)
    f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(k*x)))'
    ε_e, _ = vlasov_poisson(x, v, f₀, t)
    return t, ε_e, vmax*Δt/Δx
end

@testset "Extended verification" begin
    if get(ENV, "VASILEK_EXTENDED", "0") != "1"
        @info "extended verification skipped; set VASILEK_EXTENDED=1 to run it"
        @test true
    else
        @testset "Landau damping and dispersion" begin
            # Three wavenumbers rather than one. A single k with a single fitted
            # γ is one point on a curve, and cannot distinguish a solver that
            # reproduces the dispersion relation from one that happens to land
            # near a value at one wavenumber.
            #
            # Both roots are measured. The real frequency turns out to be the
            # sharper of the two by roughly a factor of four -- it comes from
            # counting nulls, where γ comes from fitting an amplitude that
            # numerical dissipation also acts on -- so it is held to 1% where γ
            # is held to 3%. Measured:
            #
            #   k     γ                       ω_r
            #   0.3   0.01257  (0.42%)        1.15696  (0.25%)
            #   0.4   0.06646  (0.51%)        1.28042  (0.36%)
            #   0.5   0.15558  (1.45%)        1.41372  (0.14%)
            #
            # `window` is the fitting stretch: after the initial transient, and
            # before the mode reaches the floor where recurrence and round-off
            # take over. That floor is what sets `tmax`, and it arrives *earlier
            # in units of e-foldings* the weaker the damping is -- at k = 0.3 the
            # mode has only decayed by a factor of three by t = 50, and fitting
            # past that reports γ = 0.0099, 22% low. The local rate between
            # consecutive maxima is flat at 0.0128 through t ≈ 35 and does not
            # move when Δx, Δv and Δt are all halved, so that is the estimator
            # running out of signal, not the physics changing.
            #                    k    Nx  vmax   Δt   tmax  window        alt window
            cases = ((0.5,  64,  4.0, 0.08, 70.0, (6.0, 30.0),  (8.0, 28.0)),
                     (0.4,  80,  5.0, 0.05, 60.0, (8.0, 45.0),  (10.0, 40.0)),
                     (0.3, 108,  6.0, 0.05, 60.0, (10.0, 50.0), (12.0, 45.0)))

            for (k, Nx, vmax, Δt, tmax, window, alt) in cases
                t, ε_e, courant = landau_case(k, Nx, vmax, Δt, tmax)
                γa, ωa = LANDAU_ROOTS[k]

                γ, npeaks = damping_rate(t, ε_e; tmin = window[1], tmax = window[2])
                ω, nmins  = oscillation_frequency(t, ε_e; tmin = window[1], tmax = window[2])

                println("  k = ", k, "  (Courant ", round(courant; digits = 3),
                        ", resonance at v = ", round(ωa/k; digits = 2), ")")
                println("      γ = ", rpad(round(γ; digits = 5), 8), " vs ", γa,
                        "  (", round(100*abs(γ - γa)/γa; digits = 2), "%, ",
                        npeaks, " maxima)")
                println("      ω = ", rpad(round(ω; digits = 5), 8), " vs ", ωa,
                        "  (", round(100*abs(ω - ωa)/ωa; digits = 2), "%, ",
                        nmins, " minima)")

                @test isapprox(γ, γa; rtol = 0.03)
                @test isapprox(ω, ωa; rtol = 0.01)

                # The estimator must not depend on where the window is put.
                # This is the assertion that would have caught the old fit: the
                # per-sample version moved by 2.3% when its start was nudged one
                # step, because it began on a null. Through the maxima the two
                # windows here agree to 0.43%, 0.44% and 1.45%.
                γ_alt, _ = damping_rate(t, ε_e; tmin = alt[1], tmax = alt[2])
                spread = abs(γ - γ_alt)/γ
                println("      window sensitivity: γ = ", round(γ_alt; digits = 5),
                        " on the alternate window, spread ",
                        round(100*spread; digits = 2), "%")
                @test spread < 0.03
            end
        end

        @testset "Plasma oscillations, Bohm–Gross frequency" begin
            # `docs/normalization.md` says the plasma-oscillation study verifies
            # the analytic plasma frequency. Nothing measured a frequency
            # anywhere in the repository until now -- only the energy drift
            # below, which a solver oscillating at entirely the wrong rate would
            # pass without difficulty.
            #
            # At k = 2π/100 the thermal correction is small but not negligible,
            # and that is what makes the test worth having: it separates
            # ω = √(1 + 3k²) = 1.005904 from the cold ωₚ = 1. Measured
            # ω = 1.005719, which is 0.018% from Bohm–Gross and 0.57% from the
            # cold value -- a factor of thirty, so the test fails if the thermal
            # correction is ever dropped rather than merely preferring it.
            #
            # t ≤ 200 is enough for 60 minima; the energy test below needs 3000
            # and costs forty times as much.
            x = collect(1.0:1.0:100.0)
            v = collect(-4:0.1:4)
            t = collect(0.0:0.1:200.0)
            k = 2π/100
            f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(k*x)))'
            ε_e, _ = vlasov_poisson(x, v, f₀, t)

            ω, nmins = oscillation_frequency(t, ε_e; tmin = 5.0, tmax = 195.0)
            ω_bg = sqrt(1 + 3k^2)
            println("  measured ω = ", round(ω; digits = 6),
                    "   Bohm–Gross √(1+3k²) = ", round(ω_bg; digits = 6),
                    "  (", round(100*abs(ω - ω_bg)/ω_bg; digits = 3), "%, ",
                    nmins, " minima)")
            println("  cold ωₚ = 1.0 would be off by ",
                    round(100*abs(ω - 1.0); digits = 3), "%")

            @test isapprox(ω, ω_bg; rtol = 0.002)
            # and the cold value is genuinely excluded, not merely less good
            @test abs(ω - 1.0) > 0.003
        end

        @testset "Plasma oscillations, energy conservation" begin
            x = collect(1.0:1.0:100.0)
            t = collect(0.0:0.1:3000.0)
            f(v) = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(2π*x/100)))'

            # The uniform-grid claim the notebook makes: below 0.5% at t = 3000.
            uniform = collect(-4:0.1:4)
            _, ε = vlasov_poisson(x, uniform, f(uniform), t)
            drift = (ε[end-1] - ε[1])/ε[1]
            println("  uniform grid Δε/ε = ", drift)
            @test abs(drift) < 0.005

            # The non-uniform grid was reported at "about 12%" with the limiter
            # computed globally. Per-cell-triple ξ brings it to 4.85%.
            nonuniform = vcat(collect(-4:0.2:-1.2), collect(-1:0.1:1), collect(1.2:0.2:4))
            _, ε = vlasov_poisson(x, nonuniform, f(nonuniform), t)
            drift = (ε[end-1] - ε[1])/ε[1]
            println("  non-uniform grid Δε/ε = ", drift)
            @test abs(drift) < 0.06
        end
    end
end
