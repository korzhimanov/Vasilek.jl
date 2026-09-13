# Extended verification. About a minute and a half all told -- cheap enough to
# run far more often than nightly, and CI does run it on every PR -- but still
# gated so the default test run stays instant:
#
#     VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
#
# These promote the claims the verification notebooks make in prose into
# assertions. Until they run, those claims are a 2021 HTML file.

include(joinpath(@__DIR__, "verification_harness.jl"))

# The Landau roots are **computed**, by `landau_root` in `test/dispersion.jl`,
# which the harness includes. This file used to carry three of them as typed-in
# constants; those constants now live in `test_dispersion.jl`, where the solver
# is checked against them rather than the other way round, and where they are
# the only stored numbers in a file otherwise made of identities.
#
# What this buys immediately is that `k` is no longer confined to the three
# values somebody tabulated. It also settles an old note here, which said the
# `π/(8√2)u³exp(-u²/2)` asymptotic the notebooks plot "gives 0.1145 at k = 0.5
# and is not accurate there". Three things at once, and only the conclusion was
# right: 0.1145 is that formula at `u = ω/k = 2.83`, where the notebooks
# evaluate it at `u = 1/k = 2` and get 0.3006; 0.3006 is not a `γ` but the decay
# rate of the *energy*, which is `2γ`; and as such it is 2% from the 0.30672 the
# real root gives, not 31% from anything. An asymptotic compared against the
# wrong quantity is how a 2% formula looked like a 31% error.
"`(γ, ω_r)` of the least damped Landau mode at `k`, in the order this file reads them."
function landau_rate(k)
    root = landau_root(k)
    return -imag(root), real(root)
end

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
    ε_e = vlasov_poisson(x, v, f₀, t).ε_e
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
                γa, ωa = landau_rate(k)

                γ, npeaks = damping_rate(t, ε_e; tmin = window[1], tmax = window[2])
                ω, nmins  = oscillation_frequency(t, ε_e; tmin = window[1], tmax = window[2])

                println("  k = ", k, "  (Courant ", round(courant; digits = 3),
                        ", resonance at v = ", round(ωa/k; digits = 2), ")")
                println("      γ = ", rpad(round(γ; digits = 5), 8), " vs ",
                        round(γa; digits = 5),
                        "  (", round(100*abs(γ - γa)/γa; digits = 2), "%, ",
                        npeaks, " maxima)")
                println("      ω = ", rpad(round(ω; digits = 5), 8), " vs ",
                        round(ωa; digits = 5),
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

                # And the same of ω, which rests on a counting assumption of its
                # own: the spacing is averaged between the first and last
                # minimum over `length(m)-1` intervals, so one null missed on a
                # near-tie -- or one spurious null off numerical noise --
                # rescales the answer with nothing to say so. A second window
                # changes which nulls are in the sample; a miscount does not
                # survive that, where landing near the analytic value on one
                # window could be luck.
                #
                # 1% sits in a wide gap. Below it is the estimator's own floor:
                # a null is located only to within Δt, so two windows disagree
                # by about Δt/span whatever the physics does, which is 0.4% at
                # k = 0.5. Measured here 0.10%, 0.15% and 0.05%, with a sweep
                # over further windows reaching 0.35%. Above it is the failure
                # being tested for: losing one null of ten rescales the spacing
                # by 10/9, i.e. by 11%.
                ω_alt, _ = oscillation_frequency(t, ε_e; tmin = alt[1], tmax = alt[2])
                ω_spread = abs(ω - ω_alt)/ω
                println("      window sensitivity: ω = ", round(ω_alt; digits = 5),
                        " on the alternate window, spread ",
                        round(100*ω_spread; digits = 3), "%")
                @test ω_spread < 0.01
            end
        end

        @testset "Landau damping converges under refinement" begin
            # Agreement at one resolution inside a 3% band can be luck: two
            # errors of opposite sign meeting in the middle is exactly how a
            # plausible-but-wrong solver survives a tolerance. What cannot be
            # luck is the error *shrinking* when the grid is refined.
            #
            # Δx, Δv and Δt are halved together, so the Courant number stays at
            # 0.81 and only the discretisation moves. Measured at k = 0.5:
            #
            #   Nx    Δv      Δt     γ         error    ΔL²/L² over the run
            #   32    0.2     0.16   0.16230   5.83%    -9.99e-5
            #   64    0.1     0.08   0.15558   1.45%    -9.80e-6
            #   128   0.05    0.04   0.15431   0.62%    -1.26e-6
            #   256   0.025   0.02   0.15407   0.47%    -1.62e-7
            #
            # The last level costs 14 s on its own and is left out; the three
            # below cost about two seconds together.
            #
            # The L² column is why the γ column behaves as it does, and is worth
            # asserting alongside it. `PFC` is third order, so halving the grid
            # should cut its dissipation by eight -- measured 10.2, 7.8 and 7.8.
            # The residual error in γ *is* that dissipation: the fitted rate is
            # the physical damping plus the scheme's own, which is why every
            # measurement above sits on the high side of the analytic value
            # rather than scattering about it.
            k = 0.5
            L = 2*(2π/k)
            γa = landau_rate(k)[1]
            errors = Float64[]
            dissipation = Float64[]
            for (Nx, Δv, Δt) in ((32, 0.2, 0.16), (64, 0.1, 0.08), (128, 0.05, 0.04))
                Δx = L/Nx
                x = collect(Δx:Δx:L)
                v = collect(-4:Δv:4)
                t = collect(0.0:Δt:70.0)
                f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(k*x)))'
                r = vlasov_poisson(x, v, f₀, t; invariants = true)
                γ, _ = damping_rate(t, r.ε_e; tmin = 6.0, tmax = 30.0)
                push!(errors, abs(γ - γa)/γa)
                push!(dissipation, abs(r.l2[end-1] - r.l2[1])/r.l2[1])
                println("  Nx = ", lpad(Nx, 3), "  Δv = ", rpad(Δv, 5), "  Δt = ", rpad(Δt, 4),
                        "  γ = ", rpad(round(γ; digits = 5), 7),
                        "  error ", rpad(round(100*errors[end]; digits = 2), 5), "%",
                        "  ΔL²/L² = ", round(dissipation[end]; sigdigits = 3))
            end

            # Strictly decreasing, which is the statement that distinguishes
            # convergence from coincidence.
            @test issorted(errors; rev = true)
            @test errors[end] < 0.01

            # And the dissipation falls at the scheme's order. Held loosely --
            # between 4 and 16 for a nominal 8 -- because these are three
            # coupled refinements at once and the fit window is fixed in time
            # rather than in steps, not because the measurement is noisy.
            for i in 2:length(dissipation)
                ratio = dissipation[i-1]/dissipation[i]
                println("  L² dissipation falls by ", round(ratio; digits = 1),
                        "x  (third order would give 8)")
                @test 4 < ratio < 16
            end
        end

        @testset "The invariants of the split system" begin
            # Total energy is asserted below, and was the only invariant that
            # ever was. Mass, momentum, the L² norm and the entropy are all
            # statements about the solver that the energy does not carry.
            #
            # **Measured with the cell-width sum `Σ f ΔvΔx`, not `integrate`.**
            # That is what a flux-form scheme conserves; the trapezoid halves
            # the endpoint weights, which no conservation law protects, and
            # reports 1.6e-4 of mass drift that belongs to the quadrature rather
            # than the solver. See the note on `vlasov_poisson`.
            #
            # Measured over the k = 0.5 case, 875 steps:
            #
            #   mass       2.8e-16 relative               -- round-off
            #   momentum   1.7e-15 absolute, on mass 25.1 -- round-off
            #   L²         -9.8e-6, monotone decreasing   -- numerical dissipation
            #   entropy    +7.4e-6, monotone increasing   -- the same thing
            #
            # Mass and momentum are exact conservation laws of the continuous
            # system that the discrete scheme also satisfies, so they are held
            # at round-off. L² and entropy are *not* conserved by the scheme and
            # must not be asserted as if they were: an exact Vlasov flow
            # preserves both, and the drift is the numerical dissipation the
            # refinement test above measures. What can be asserted is the
            # **direction** -- a dissipative scheme can only lose L² and gain
            # entropy, and a step of either the wrong way is a bug.
            k = 0.5
            L = 2*(2π/k)
            Δx = L/64
            x = collect(Δx:Δx:L)
            v = collect(-4:0.1:4)
            t = collect(0.0:0.08:70.0)
            f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(k*x)))'
            r = vlasov_poisson(x, v, f₀, t; invariants = true)

            # `[end]` is padded from `[end-1]` by the driver, so the histories
            # are read one short throughout.
            mass, mom = r.mass[1:end-1], r.momentum[1:end-1]
            l2, S = r.l2[1:end-1], r.entropy[1:end-1]

            mass_drift = maximum(abs, mass .- mass[1])/mass[1]
            println("  mass     drift ", mass_drift, " relative")
            @test mass_drift < 1e-13

            # Momentum starts at zero by symmetry, so there is no scale to take
            # a ratio against; it is held in absolute terms against the mass.
            println("  momentum max|p| ", maximum(abs, mom), "  on mass ", round(mass[1]; digits = 3))
            @test maximum(abs, mom) < 1e-12

            println("  L²       ", (l2[end] - l2[1])/l2[1], " relative, monotone: ",
                    all(diff(l2) .≤ 1e-15))
            @test l2[end] < l2[1]
            @test all(diff(l2) .≤ 1e-15)          # dissipative at every step

            println("  entropy  ", (S[end] - S[1])/S[1], " relative, monotone: ",
                    all(diff(S) .≥ -1e-15))
            @test S[end] > S[1]
            @test all(diff(S) .≥ -1e-15)          # and never un-mixes
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
            ε_e = vlasov_poisson(x, v, f₀, t).ε_e

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

            # The discrimination above is 0.57% wide, so it is worth knowing
            # that ω is not an artefact of its window. Same argument as in the
            # Landau loop, and cheap here because no second run is needed --
            # only a second reading of the same ε_e. The long window makes this
            # the sharpest of the four: 60 minima over 190 time units puts the
            # Δt/span floor at 0.05%, and a sweep of windows moves ω by at most
            # 0.028%, where losing one null of sixty would move it by 1.7%.
            ω_alt, nalt = oscillation_frequency(t, ε_e; tmin = 20.0, tmax = 180.0)
            ω_spread = abs(ω - ω_alt)/ω
            println("  window sensitivity: ω = ", round(ω_alt; digits = 6),
                    " on t ∈ [20, 180] (", nalt, " minima), spread ",
                    round(100*ω_spread; digits = 4), "%")
            @test ω_spread < 0.002
        end

        @testset "Plasma oscillations, energy conservation" begin
            x = collect(1.0:1.0:100.0)
            t = collect(0.0:0.1:3000.0)
            f(v) = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(2π*x/100)))'

            # The uniform-grid claim the notebook makes: below 0.5% at t = 3000.
            uniform = collect(-4:0.1:4)
            ε = vlasov_poisson(x, uniform, f(uniform), t).ε
            drift = (ε[end-1] - ε[1])/ε[1]
            println("  uniform grid Δε/ε = ", drift)
            @test abs(drift) < 0.005

            # The non-uniform grid was reported at "about 12%" with the limiter
            # computed globally. Per-cell-triple ξ brings it to 4.85%.
            nonuniform = vcat(collect(-4:0.2:-1.2), collect(-1:0.1:1), collect(1.2:0.2:4))
            ε = vlasov_poisson(x, nonuniform, f(nonuniform), t).ε
            drift = (ε[end-1] - ε[1])/ε[1]
            println("  non-uniform grid Δε/ε = ", drift)
            @test abs(drift) < 0.06
        end

        @testset "Two-stream instability" begin
            # The first *unstable* case in the suite. Everything else here is a
            # damped or neutral mode, and a growth rate catches a class of error
            # that damping cannot: a sign flip in the field push turns damping
            # into growth and growth into damping, so a suite made only of
            # damped cases is half-blind to it.
            #
            # **The growth rate is derived here, not quoted.** For two cold
            # beams of density 1/2 at ±v₀ the electrostatic dispersion relation
            # is
            #
            #     1 = ½/(ω − kv₀)² + ½/(ω + kv₀)²
            #
            # which, with `a = kv₀` and `u = ω²`, is a quadratic:
            #
            #     u² − (2a² + 1)u + (a⁴ − a²) = 0
            #     u± = [(2a² + 1) ± √(8a² + 1)]/2
            #
            # `u₋ < 0` exactly when `a < 1`, and then `γ = √(−u₋)`. So this
            # needs no plasma dispersion function, no numerical root, and no
            # tabulated constant of the kind the Landau cases have to carry --
            # which is what made the case worth waiting for rather than
            # hard-coding a number.
            #
            # **The beams this runs are warm, and the comparison is now against
            # warm theory.** `γ_cold` stays as the `vt → 0` limit -- checked as
            # such in `test_dispersion.jl` -- but it is not what the simulation
            # is held to any more. `two_stream_warm` solves the same relation
            # with Maxwellian beams instead of delta functions, and the
            # difference is not cosmetic: at the `vt = 0.3` these runs use the
            # cold form is off by up to 7.44% where the warm root is off by
            # 2.01%, which is the whole tolerance budget. The one place it
            # changes a *conclusion* rather than a number is `a = 0.8`; see
            # below.

            @testset "the closed form solves the dispersion relation" begin
                # Cheap, and it is what lets the rest of this testset be read as
                # a measurement rather than a comparison against a constant
                # somebody typed. Measured residual ≤ 1.4e-14.
                worst = 0.0
                for a in (0.2, 0.4, 0.6, 0.8, 0.95)
                    worst = max(worst, abs(two_stream_residual(im*γ_cold(a), a)))
                end
                println("  closed form vs dispersion relation: worst |residual| = ", worst)
                @test worst < 1e-12

                # and it reproduces the two textbook numbers on its own: the
                # fastest-growing wavenumber sits at a = √(3/8) and the peak
                # rate at 1/(2√2).
                @test isapprox(γ_cold(sqrt(3/8)), 1/(2*sqrt(2)); rtol = 1e-12)
                grid = 0.001:0.001:0.999
                @test isapprox(grid[argmax(γ_cold.(grid))], sqrt(3/8); atol = 2e-3)
                # stable above the boundary, exactly
                @test γ_cold(1.0) == 0.0
                @test γ_cold(1.5) == 0.0
            end


            @testset "the growth rate follows the dispersion relation" begin
                # Three wavenumbers, not one -- and `γ(a)` is **non-monotone**,
                # rising to a peak near a = √(3/8) ≈ 0.612 (exactly there for
                # the cold branch, close to it for the warm one) and falling
                # again, so reproducing all three is a statement about the
                # branch rather than about one point. A solver that merely
                # amplified whatever it was given could not put the maximum in
                # the right place.
                #
                # Measured, with the beams at vt = 0.3 and the fit taken over
                # ε_e rising from 100x its initial value to 5.0:
                #
                #   a = kv₀   γ measured   γ warm     error    γ cold     error
                #   0.4       0.30244      0.30362   -0.39%    0.30819   -1.87%
                #   0.6       0.34229      0.34909   -1.95%    0.35339   -3.14%
                #   0.8       0.31232      0.31201   +0.10%    0.31134   +0.31%
                #
                # Held to 3% against the warm root, about 1.5x the worst. The
                # cold column is printed alongside because the two disagree by
                # more than the tolerance and it is worth seeing which one the
                # solver follows.
                #
                # **`a = 0.8` is not an overshoot.** It reads +0.31% above the
                # cold value, and an earlier version of this comment attributed
                # that to the beat ripple on the grounds that finite temperature
                # cannot make a beam grow faster than a cold one. That premise
                # is false near the band edge: the warm rate crosses above the
                # cold one between a = 0.75 and a = 0.8, and beyond a = 1 the
                # cold branch is identically zero while warm beams are still
                # unstable -- which this very testset measures below. Against
                # the warm root the same measurement is +0.10%, and nothing
                # needs explaining away.
                measured = Float64[]
                for a in (0.4, 0.6, 0.8)
                    t, ε_e = two_stream(a)
                    γ, t0, t1 = growth_rate(t, ε_e; lo = 100*ε_e[1], hi = 5.0)
                    γw = two_stream_warm(a)
                    push!(measured, γ)
                    println("  a = kv₀ = ", a, "  γ = ", round(γ; digits = 5),
                            " vs warm ", round(γw; digits = 5),
                            " (", round(100*(γ - γw)/γw; digits = 2), "%)",
                            ", cold ", round(γ_cold(a); digits = 5),
                            " (", round(100*(γ - γ_cold(a))/γ_cold(a); digits = 2),
                            "%), fitted over t ∈ [", round(t0; digits = 2), ", ",
                            round(t1; digits = 2), "]")
                    @test isapprox(γ, γw; rtol = 0.03)
                end

                # The shape, independent of the individual tolerances: the
                # middle wavenumber is the fastest-growing one.
                @test measured[2] > measured[1]
                @test measured[2] > measured[3]

                # Temperature at fixed wavenumber, where the effect is
                # unambiguous: at a = 0.6 the warm branch falls monotonically
                # with vt (0.35337 at 0.02, 0.34909 at 0.3, 0.33381 at 0.6 --
                # `test_dispersion.jl` asserts the monotonicity), and the runs
                # follow it. Measured γ = 0.32710 at vt = 0.6 against 0.34229 at
                # vt = 0.3, a drop of 4.4% where theory predicts 4.4%.
                #
                # Both halves are asserted now. The relative one -- a colder
                # beam grows faster -- is the statement that survives any error
                # common to the two runs. The absolute one is possible only
                # because the target moved with the temperature: against the
                # cold value this run is 7.44% out, which is why the old comment
                # could compare the two runs with each other but not either with
                # theory.
                t_wide, ε_wide = two_stream(0.6; vt = 0.6)
                γ_wide, _, _ = growth_rate(t_wide, ε_wide; lo = 100*ε_wide[1], hi = 5.0)
                γw_wide = two_stream_warm(0.6; vt = 0.6)
                println("  vt = 0.6 gives γ = ", round(γ_wide; digits = 5),
                        " vs warm ", round(γw_wide; digits = 5),
                        " (", round(100*(γ_wide - γw_wide)/γw_wide; digits = 2), "%)",
                        ", against ", round(measured[2]; digits = 5), " at vt = 0.3",
                        "  (cold, for both, ", round(γ_cold(0.6); digits = 5), ")")
                @test γ_wide < measured[2]
                @test isapprox(γ_wide, γw_wide; rtol = 0.03)
            end

            @testset "and stops at the stability boundary" begin
                # Both forms of the theory make the cases below stable, and they
                # do it at different places: `γ_cold` is exactly zero for
                # `kv₀ ≥ 1`, while the warm band edge at vt = 0.3 sits between
                # a = 1.0 and a = 1.05 -- `two_stream_warm` returns 0.09823 at
                # 1.0 and 0.0 at 1.05. a = 1.2 and 1.6 are therefore stable by
                # both, which is what makes them the qualitative cases: no
                # tolerance can launder a failure. A solver with the field sign
                # reversed, or one amplifying grid noise, grows here.
                @test two_stream_warm(1.2) == 0.0
                @test two_stream_warm(1.6) == 0.0

                # Measured over t ≤ 26, as a ratio of peak ε_e to initial:
                # a = 1.2 gives 1.00 (4.84e-5 decaying to 3.36e-5) and a = 1.6
                # gives 1.00 (2.02e-5 to 7.02e-6), against 5.9e4 at a = 0.6.
                for a in (1.2, 1.6)
                    t, ε_e = two_stream(a; tmax = 26.0)
                    ratio = maximum(ε_e)/ε_e[1]
                    println("  a = kv₀ = ", a, " (stable): peak/initial ε_e = ",
                            round(ratio; digits = 3), ", final/initial = ",
                            round(ε_e[end]/ε_e[1]; digits = 3))
                    @test ratio < 1.5
                    @test ε_e[end] < ε_e[1]        # Landau-damped, not merely flat
                end

                # `a = 1.0` is the cold boundary itself, where the cold form
                # predicts exactly nothing and the warm one predicts 0.09823.
                # It used to be asserted as *growth of some kind*, a factor of
                # 4.91 over t ≤ 26, because a rate could not be compared against
                # anything: the only theory in the file said zero. With the warm
                # root it becomes the sharpest case in this testset -- a number
                # against a number, in a regime where the two theories differ by
                # infinity rather than by percent.
                #
                # Measured γ = 0.09510 against 0.09823, which is 3.19%, fitted
                # over t ∈ [45.6, 78.3]. The band is the same `[100ε₀, 5.0]`
                # every other case uses; `tmax = 80` is what it takes to reach
                # 5.0 at a tenth of the growth rate, and the run goes non-finite
                # at t = 86.3, so the margin is six time units rather than the
                # four steps the `a = 0.6` case gets at `tmax = 24`.
                #
                # Held to 8% rather than the 3% above, and the reason is in
                # `growth_rate`: the beat between the growing root and the
                # oscillating pair decays as exp(-γt) relative to the mode, so
                # it is worst where γ is smallest, and γ here is a third of the
                # branch maximum. Measured over `hi` ∈ {1, 2, 3, 5} the fit
                # moves from -6.08% to -3.19%; 8% covers that spread with room.
                t, ε_e = two_stream(1.0; tmax = 80.0)
                γ, t0, t1 = growth_rate(t, ε_e; lo = 100*ε_e[1], hi = 5.0)
                γw = two_stream_warm(1.0)
                println("  a = kv₀ = 1.0 (the cold boundary): γ = ", round(γ; digits = 5),
                        " vs warm ", round(γw; digits = 5),
                        " (", round(100*(γ - γw)/γw; digits = 2), "%), cold says 0",
                        ", fitted over t ∈ [", round(t0; digits = 1), ", ",
                        round(t1; digits = 1), "]")
                @test isapprox(γ, γw; rtol = 0.08)
                @test all(isfinite, ε_e)
            end
        end

        @testset "Laser wakefield: the laser drives a plasma wave" begin
            # This testset used to assert that the study "runs and stays
            # bounded", because that was all there was to assert: `wakefield`
            # had no ponderomotive coupling, so the laser never entered the
            # longitudinal push and `ex` was the slab edges relaxing. The wake
            # and energy numbers came out bit-identical whether the transverse
            # current was right, wrong by a factor of `Δt`, or wrong by
            # thirty-two orders of magnitude -- a test that cannot see the laser
            # is not a test of a laser wakefield.
            #
            # What follows measures the wake against linear wakefield theory
            # instead. The claims are in four groups, and they fail
            # independently: that a plasma wave is there at all (frequency,
            # wavelength), that the *laser* made it (phase locking, causality,
            # the `a₀²` law, the unlit control), that it is the right size
            # (amplitude and pointwise profile against `linear_wake`), and that
            # the solver stayed sane while doing it.
            r = wakefield()
            n₀, T = r.plasma_density, r.plasma_temperature
            ωₚ = sqrt(n₀)
            x, t = r.x, r.t

            # The slab interior at the final time: inside `[0, 10·2π]`, clear of
            # the sheaths the two edges carry. The control run below puts those
            # at 2.2e-3 within two units of each edge and 2.8e-4 in here, which
            # is the 3% floor under every comparison in this testset.
            window = findall(i -> 8.0 ≤ x[i] ≤ 55.0, eachindex(x))

            # ---- the driver, measured *and* predicted
            #
            # It used to be only measured, on the grounds that a one-cycle
            # driver has no single group velocity and there was nothing to
            # predict it with. `vg_pulse` is that closed form -- the group
            # velocity of the *discrete* dispersion relation, averaged over the
            # pulse's own spectrum -- and it turns the driver from an input of
            # this testset into one more thing under test.
            #
            # Measured 0.9261 against 0.9329, which is 0.73%. Held to 2%, the
            # same tolerance as the phase-locking claim below, because the
            # residue is of the same kind: `pulse_velocity` fits the position of
            # a peak that is quantised to Δx and slips by half a carrier
            # wavelength whenever the envelope's maximum crosses a crest.
            Δx_run, Δt_run = x[2] - x[1], t[2] - t[1]
            v, nv = pulse_velocity(t, x, r.Φ; lo = 5.0, hi = 55.0)
            v_theory = vg_pulse(n₀, Δx_run, Δt_run; duration = 2π)
            println("  driver: measured ", round(v; digits = 4), " vs vg_pulse ",
                    round(v_theory; digits = 4),
                    " (", round(100*(v - v_theory)/v_theory; digits = 2), "%)",
                    ", continuum √(1-n) = ", round(sqrt(1 - n₀); digits = 4))
            @test isapprox(v, v_theory; rtol = 0.02)

            # ---- a plasma wave is there
            λ, nλ = wave_period(x, r.ex[end, :]; lo = 8.0, hi = 55.0)
            k = 2π/λ
            i₂₅ = argmin(abs.(x .- 25.0))
            # Start the time series two pulse durations after the pulse passed
            # this point, which the `Φ` history dates rather than the defaults.
            t₀ = t[argmax(view(r.Φ, :, i₂₅))] + 2*2π
            period, nper = wave_period(t, r.ex[:, i₂₅]; lo = t₀, hi = t[end])
            ω = 2π/period

            # The wake oscillates at the plasma frequency, Bohm--Gross included:
            # measured 0.32121 against `√(ωₚ² + 3Tk²)` = 0.32195, which is
            # 0.23%, and against the cold ωₚ = 0.31623, which is 1.58%. The warm
            # value is the better fit and is the one asserted, but the two are
            # only 1.8% apart and the measurement moves by 0.5% with the probe
            # point, so -- unlike the standing oscillation above, which resolves
            # the correction at 0.018% against 0.57% -- this run is not entitled
            # to claim the cold value is excluded. It is asserted against the
            # right theory; it does not discriminate between the two.
            #
            # This is also the one number here that does *not* move with
            # resolution: 0.32116 at half the cells against 0.32121 here. The
            # frequency is set by the plasma, which both grids resolve; it is the
            # wavelength and the driver's speed that the laser's grid controls.
            ω_bg = sqrt(ωₚ^2 + 3*T*k^2)
            println("  pulse v = ", round(v; digits = 4), " ($nv samples)",
                    "   λ = ", round(λ; digits = 3), " ($nλ zeros)",
                    " against ", round(wake_wavelength(v, n₀, T); digits = 3),
                    "   ω = ", round(ω; digits = 5), " ($nper zeros)",
                    " against Bohm-Gross ", round(ω_bg; digits = 5),
                    " (cold ωₚ = ", round(ωₚ; digits = 5), ")")
            @test isapprox(ω, ω_bg; rtol = 0.02)

            # The wavelength is the driver's own: a wave that keeps station with
            # something moving at `v` has `ω(k) = kv`, which against Bohm--Gross
            # gives `2π√(v² - 3T)/ωₚ`. Measured 18.010 against 18.076 from the
            # measured driver, 0.37%, and stable to 0.4% across the windows
            # [8,42], [8,55] and [5,58].
            @test isapprox(λ, wake_wavelength(v, n₀, T); rtol = 0.03)

            # And against the *predicted* driver, which closes the loop: no
            # quantity measured in this run enters the right-hand side, so the
            # wavelength is now a prediction from the grid parameters and the
            # plasma alone. Measured 18.010 against 18.214, 1.12% -- larger than
            # the 0.37% above, as it must be, since it carries the driver's own
            # 0.73% as well.
            @test isapprox(λ, wake_wavelength(v_theory, n₀, T); rtol = 0.03)

            # ---- the laser made it
            #
            # Phase locking is the statement here that cannot come from the
            # slab: the wave's own phase velocity `ω/k`, from two independent
            # measurements along two different axes, is the speed of the pulse.
            # Measured 0.9207 against 0.9261, 0.58%.
            #
            # Both sit below the continuum `√(1-n)` = 0.9487, and the gap is the
            # grid rather than the physics -- `vg_pulse` predicts 0.9329 here.
            # At half the resolution the same three numbers read 0.8925, 0.8858
            # and 0.8993: the phase locking holds just as well while the thing
            # being locked to is six percent slow. That is why the driver is now
            # asserted against a closed form as well, and why the study runs at
            # twenty cells per wavelength.
            @test isapprox(ω/k, v; rtol = 0.02)

            # Behind the pulse and not ahead of it. The margin is two pulse
            # durations, where the drive is down to 1e-4 of its peak -- at one
            # margin the window is still inside the pulse, the Gaussian being as
            # wide as the wake is long. Measured 9.4, and the field that is
            # ahead is not all precursor: the Poisson solve is instantaneous, so
            # the charge bunches behind do reach forward, and the unlit control
            # leaves a wake of 3.4e-4 without any laser at all.
            kmid = findfirst(k -> x[argmax(view(r.Φ, k, :))] ≥ 30.0 &&
                                  r.Φ[k, argmax(view(r.Φ, k, :))] > 0.5*maximum(r.Φ),
                             1:length(t))
            # Guarded the way the harness guards its own lookups: unguarded,
            # a run whose pulse never reaches x = 30 at half its peak -- a
            # shorter `total_time`, a moved slab -- fails as a `MethodError`
            # inside `view(Φ, nothing, :)` rather than saying what is missing.
            kmid === nothing &&
                error("the pulse never reached x = 30 above half its peak Φ; " *
                      "there is no mid-slab snapshot to compare behind against ahead")
            ipk = argmax(view(r.Φ, kmid, :))
            behind = findall(i -> 5.0 ≤ x[i] ≤ x[ipk] - 2*2π, eachindex(x))
            ahead = findall(i -> x[ipk] + 2*2π ≤ x[i] ≤ 58.0, eachindex(x))
            (isempty(behind) || isempty(ahead)) &&
                error("the pulse at x = $(x[ipk]) leaves no room for a two-duration " *
                      "margin on both sides inside [5, 58]")
            println("  at t = ", round(t[kmid]; digits = 1), " the pulse is at x = ",
                    round(x[ipk]; digits = 1), ": |ex| behind = ",
                    round(maximum(abs, r.ex[kmid, behind]); digits = 6), ", ahead = ",
                    round(maximum(abs, r.ex[kmid, ahead]); digits = 6))
            @test maximum(abs, r.ex[kmid, behind]) >
                  4*maximum(abs, r.ex[kmid, ahead])

            # ---- the right size, against linear theory
            #
            # `linear_wake` integrates the driven plasma oscillator on the `Φ`
            # this run recorded. It shares no code with the wake it is compared
            # against: `Φ` is the `FDTD1D` side, `ex` is the Vlasov push and the
            # Poisson solve, and the ODE is the only thing claiming they agree.
            ref = linear_wake(x, t, r.Φ, r.nᵢ; temperature = T)
            amp = maximum(abs, r.ex[end, window])
            amp_ref = maximum(abs, ref[window])
            rms = sqrt(sum((r.ex[end, window] .- ref[window]).^2)/sum(ref[window].^2))
            println("  wake amplitude = ", round(amp; digits = 6), " against ",
                    round(amp_ref; digits = 6), " from linear theory (",
                    round(100*(amp/amp_ref - 1); digits = 2), "%), pointwise rms ",
                    round(rms; digits = 4))

            @test isapprox(amp, amp_ref; rtol = 0.10)   # measured 6.4% under

            # Pointwise, not just in amplitude: the profiles agree to 6.4% of
            # the theory's own rms over the window, which is phase as well as
            # size. Dropping the `3T` term from `linear_wake` -- the thermal
            # correction alone, worth 1.7% on the frequency -- takes this to
            # 0.32 and fails, so the tolerance is not loose enough to pass a
            # reference with the physics wrong. It shows up here rather than in
            # the amplitude, which that same change moves by 4%, well inside the
            # 10% above.
            @test rms < 0.15

            # ---- and it is the laser's, at the laser's own scaling
            #
            # The wake of a ponderomotive drive goes as the intensity, so
            # halving the amplitude quarters it. Measured ratio 4.13 against 4,
            # which is 3.3% -- and it is the one claim here that a compensating
            # error in the drive and in the response cannot fake, since it holds
            # the plasma fixed and moves only the laser.
            #
            # **The pair is 0.3 and 0.15 because further down the ladder the
            # measurement stops being about the laser.** Continuing the sweep:
            #
            #   cells/λ₀   a₀ = 0.3   0.15      0.075     0.3/0.15   0.15/0.075
            #   20         0.012366   0.002993  0.000788  4.13       3.80
            #   10         0.008900   0.002235  0.000671  3.98       3.33
            #
            # The unlit control leaves a wake of 3.4e-4 -- the slab edges
            # relaxing -- which is 3% of the a₀ = 0.3 amplitude and 43% of the
            # a₀ = 0.075 one, so the bottom rung is measuring the sheaths as much
            # as the laser and its ratio drifts accordingly. The few percent left
            # at the top of the ladder is not attributed here; it is smaller than
            # the spread between resolutions, which is the honest bound on it.
            half = wakefield(laser_amplitude = 0.15)
            ratio = amp/maximum(abs, half.ex[end, window])
            println("  amplitude ratio at half the laser = ", round(ratio; digits = 4))
            @test isapprox(ratio, 4.0; rtol = 0.08)

            # The control: the same run with the laser off. This is what the old
            # testset was measuring without knowing it -- the slab edges
            # relaxing -- and the wake is 36 times it. Before the ponderomotive
            # term the ratio here would have been 1.
            dark = wakefield(laser_amplitude = 0.0)
            println("  wake with the laser off = ",
                    round(maximum(abs, dark.ex[end, window]); digits = 6),
                    ", lit/unlit = ",
                    round(amp/maximum(abs, dark.ex[end, window]); digits = 1))
            @test amp > 10*maximum(abs, dark.ex[end, window])

            # ---- and the solver stayed sane
            @test all(isfinite, r.ey)
            @test all(isfinite, r.ex)
            @test all(isfinite, r.n)

            # `PFC` is positivity preserving and the density is its momentum
            # integral, so this must hold exactly. Measured minimum: 0.0.
            @test minimum(r.n) ≥ 0.0

            # `ε` is the longitudinal energy: `∫∫f p² + ∫e²`, kinetic *and*
            # electrostatic. Both halves matter to what a rise means. Because
            # the field term is in, a plasma oscillation trading kinetic energy
            # for field energy leaves `ε` alone, so it is a genuine longitudinal
            # invariant and not a quantity that sloshes on its own. Because the
            # transverse motion and the transverse field are out, a laser doing
            # work on the plasma pushes energy across that boundary and `ε` is
            # *supposed* to rise -- 13.7% here against the 1.2% it drifted when
            # nothing was coupled, and against 8.4% at half this resolution,
            # where the wake it drives is 39% smaller. It is a bound against
            # divergence, not a conservation claim, and it is the one number in
            # this testset that predicts nothing.
            println("  Δε/ε = ", round((r.ε[end] - r.ε[1])/r.ε[1]; digits = 4),
                    "   peak |p⊥| = ", round(sqrt(2*maximum(r.Φ)); digits = 4),
                    "   min n = ", minimum(r.n))
            @test 0 < (r.ε[end] - r.ε[1])/r.ε[1] < 0.2

            # ---- and the driver's error is the grid's, which refining shows
            #
            # Everything above is one resolution, and one resolution cannot
            # distinguish "the driver moves at 0.93" from "the driver moves at
            # 0.95 and this grid is 2% slow". Refining is what separates them,
            # and it is the reason this study runs at twenty cells per
            # wavelength rather than the ten it used to.
            #
            #   cells/λ₀   measured v   vg_pulse   error    λ        Δε/ε
            #   10         0.8858       0.8993     -1.51%   17.460   0.084
            #   20         0.9261       0.9329     -0.73%   18.010   0.137
            #
            # Two statements, and they are different. The first is that each
            # measurement matches the closed form *for its own grid*: the
            # discrepancy is the estimator, not the physics. The second is that
            # the sequence moves toward the continuum √(1-n) = 0.9487 rather
            # than toward wherever the coarse grid happened to sit -- which is
            # what says the closed form is the right one and not a curve fitted
            # to one run.
            #
            # The coarse run costs 4 s. Nothing else in this testset is repeated
            # at it: the wake's *frequency* barely moves (0.32116 against
            # 0.32121), and the comparison against `linear_wake` cannot see the
            # grid at all, because that reference is driven by the `Φ` of the run
            # it is checking -- both sides shift together. That is a strength for
            # what it does test and a blind spot for this, and it is why the
            # driver needs a closed form of its own.
            coarse = wakefield(Δx = 0.1*2π)
            Δx_c, Δt_c = coarse.x[2] - coarse.x[1], coarse.t[2] - coarse.t[1]
            v_c, _ = pulse_velocity(coarse.t, coarse.x, coarse.Φ; lo = 5.0, hi = 55.0)
            v_c_theory = vg_pulse(n₀, Δx_c, Δt_c; duration = 2π)
            λ_c, _ = wave_period(coarse.x, coarse.ex[end, :]; lo = 8.0, hi = 55.0)
            println("  at ", round(Int, 2π/Δx_c), " cells/λ₀: v = ", round(v_c; digits = 4),
                    " vs vg_pulse ", round(v_c_theory; digits = 4),
                    " (", round(100*(v_c - v_c_theory)/v_c_theory; digits = 2), "%)",
                    ", λ = ", round(λ_c; digits = 3),
                    "   against ", round(Int, 2π/Δx_run), " cells: ",
                    round(v; digits = 4), ", ", round(λ; digits = 3))
            @test isapprox(v_c, v_c_theory; rtol = 0.02)
            @test abs(v - sqrt(1 - n₀)) < abs(v_c - sqrt(1 - n₀))
        end
    end
end
