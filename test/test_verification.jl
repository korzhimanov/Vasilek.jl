# Extended verification. About a minute and a half all told -- cheap enough to
# run far more often than nightly, and CI does run it on every PR -- but still
# gated so the default test run stays instant:
#
#     VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
#
# These promote the claims the verification notebooks make in prose into
# assertions. Until they run, those claims are a 2021 HTML file.

@isdefined(vlasov_poisson) || include(joinpath(@__DIR__, "verification_harness.jl"))

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
    trapping_phase(α, γ, T)

Bounce phase a resonant particle accumulates before the mode damps away,
`∫₀ᵀ ω_B dt` with `ω_B = √(kE) = √α` falling as `exp(-γt)` with the field:

    Φ_B = (√α/γ)·(1 - exp(-γT))

**This is the number that says whether a "linear" Landau run is linear.** The
usual statement is O'Neil's `γ ≫ ω_B`, and by that measure none of the cases
below qualify -- but `ω_B` is not constant, and a mode that damps before the
first bounce never traps anything however small `γ` is. The phase is what
combines the two, and it is calibrated against the trapping testset below --
**in its own units**, which matters because that testset reports the undamped
phase `√α·t₀` (7.1 to 8.0 at the arrest) rather than this one. Evaluated at the
same arrest times, `Φ_B` is 4.3 to 5.7; at the point where the local rate first
leaves the analytic one (`t ≈ 35` at `α = 0.01`) it is 2.8. The guard the linear
cases use, `Φ_B < 2`, sits below both, and the testset asserts it against the
arrest.

Measured for the three cases here, at `α = 1e-3`: 0.21, 0.45 and 1.70. For the
`α = 0.01` this suite used at `k = 0.3`, over its old window: 3.7 -- which is
where the damping had started to stop, and the reason that case read 0.42% *low*
where every other measurement in the suite reads high.
"""
trapping_phase(α, γ, T) = sqrt(α)/γ*(1 - exp(-γ*T))

"""
    landau_case(k, Nx, vmax, Δt, tmax; α = 1e-3)

A single-mode Landau run: two wavelengths of `k` across the box, a Maxwellian
perturbed by `α`, returning `(t, ε_e, courant, α)`. The amplitude comes back with
the run so that the linearity guard reads it off the result rather than
restating it -- restated, the guard was a check on a constant, and passed at the
`α = 1e-2` it exists to catch.

**`α = 1e-3`, not the 1e-2 this suite started with.** The mode is meant to be
linear, and at 1% it is not: see [`trapping_phase`](@ref) and the trapping
testset. Dropping the amplitude costs nothing in signal -- the fit runs through
the maxima of a quantity that spans decades either way -- and it moves the
`k = 0.3` measurement from 0.42% below the analytic rate to 0.71% above it,
which is where the other two sit and where numerical dissipation puts them.

**The grid is not free.** The velocity window has to contain the resonance at
`v = ω_r/k` -- 2.83, 3.21 and 3.87 for the three cases below -- because that is
where the damping comes from; and the window then fixes the time step, since the
fastest row runs at `max|v|·Δt/Δx` and `PFCNonUniform` is a finite-volume scheme
with a Courant limit of 1. Widening the window to reach a resonance therefore
costs a smaller `Δt`, not just more velocity points.

The number returned and printed as the Courant figure is `max|v|·Δt/Δx`, which
runs between 0.64 and 0.82 here. The *sweeps* run at half of it: Strang takes
two x half-steps per step, so the displacement per call is `vΔt/2` and the real
figure is 0.32 to 0.41. The conservative number is the one to reason with when
widening a velocity window, which is why it is the one shown.
"""
function landau_case(k, Nx, vmax, Δt, tmax; α = 1e-3)
    L = 2*(2π/k)
    Δx = L/Nx
    # Nx points by construction. `Δx:Δx:L` rounds its length out of the
    # endpoints and can come up a point short; see `two_stream`.
    x = collect(range(Δx; step = Δx, length = Nx))
    v = collect(-vmax:0.1:vmax)
    t = collect(0.0:Δt:tmax)
    f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + α*cos(k*x)))'
    ε_e = vlasov_poisson(x, v, f₀, t).ε_e
    return t, ε_e, vmax*Δt/Δx, α
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
            # sharper of the two by a factor of four to nine -- it comes from
            # counting nulls, where γ comes from fitting an amplitude that
            # numerical dissipation also acts on -- so it is held to 1% where γ
            # is held to 3%. Measured:
            #
            #   k     γ                       ω_r
            #   0.3   0.01271  (0.71%)        1.15895  (0.08%)
            #   0.4   0.06688  (1.14%)        1.28228  (0.22%)
            #   0.5   0.15481  (0.94%)        1.41939  (0.26%)
            #
            # All three sit *above* the analytic rate now, and that is the point
            # of the amplitude change: the residue is numerical dissipation,
            # which can only add damping. At `α = 1e-2` the k = 0.3 case read
            # 0.42% below instead -- trapping pulling one way while dissipation
            # pulled the other, and the agreement was the two cancelling. See
            # [`trapping_phase`](@ref) and the trapping testset.
            #
            # `window` is the fitting stretch: after the initial transient, and
            # before the mode reaches the floor where recurrence and round-off
            # take over. At k = 0.3 the damping is slow enough that the window
            # can run to t = 90 -- thirty maxima -- which at α = 1e-2 would have
            # been deep into the trapped regime.
            #                    k    Nx  vmax   Δt    tmax   window        alt window
            cases = ((0.5,  64,  4.0, 0.08,  70.0, (6.0, 30.0),  (8.0, 28.0)),
                     (0.4,  80,  5.0, 0.05,  60.0, (8.0, 45.0),  (10.0, 40.0)),
                     (0.3, 108,  6.0, 0.05, 100.0, (10.0, 90.0), (12.0, 80.0)))

            for (k, Nx, vmax, Δt, tmax, window, alt) in cases
                t, ε_e, courant, α = landau_case(k, Nx, vmax, Δt, tmax)
                γa, ωa = landau_rate(k)

                # The run has to be linear for the analytic rate to be the right
                # target, and that is a property of the amplitude and the window
                # together rather than of either alone. Measured 0.21, 0.45 and
                # 1.70 for the three cases; in these units the damping departs at
                # about 2.8 and stops at 4.3 to 5.7. `α` is the run's own: at the
                # 1e-2 this suite used to run, the k = 0.3 case reads 5.38 here.
                Φ_B = trapping_phase(α, γa, window[2])
                @test Φ_B < 2.0

                γ, npeaks = damping_rate(t, ε_e; tmin = window[1], tmax = window[2])
                ω, nmins  = oscillation_frequency(t, ε_e; tmin = window[1], tmax = window[2])

                println("  k = ", k, "  (Courant ", round(courant; digits = 3),
                        ", resonance at v = ", round(ωa/k; digits = 2),
                        ", bounce phase ", round(Φ_B; digits = 2), ")")
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
                # windows here agree to 0.04%, 0.05% and 0.18% -- tighter than
                # the 0.43%, 0.44% and 1.45% of the 1% runs, the k = 0.3 case by
                # a factor of eight, because what moved that one between windows
                # was trapping rather than the estimator.
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
                # k = 0.5. Measured here 0.300%, 0.000% and 0.066%. Above it is
                # the failure being tested for: losing one null of ten rescales
                # the spacing by 10/9, i.e. by 11%.
                ω_alt, _ = oscillation_frequency(t, ε_e; tmin = alt[1], tmax = alt[2])
                ω_spread = abs(ω - ω_alt)/ω
                println("      window sensitivity: ω = ", round(ω_alt; digits = 5),
                        " on the alternate window, spread ",
                        round(100*ω_spread; digits = 3), "%")
                @test ω_spread < 0.01
            end
        end

        @testset "Trapping stops the damping, on the bounce time" begin
            # The linear cases above are linear because their amplitude was
            # chosen to make them so. This is the same solver at amplitudes
            # where it is not, and it is a physics test rather than a caveat:
            # the arrest of Landau damping by trapped particles is O'Neil's
            # result, and reproducing *when* it happens is a statement about the
            # nonlinear term that no damped-mode test can make.
            #
            # The local rate between maxima two apart, at k = 0.3, α = 1e-2:
            #
            #   t      8     19    30    41    52    63    73    84    95
            #   γ      .0127 .0130 .0128 .0119 .0092 .0048 .0000 -.0036 -.0055
            #
            # It does not merely stop -- it goes negative, which is the field
            # growing again as the trapped population sloshes. At α = 1e-3 the
            # same column is flat at 0.0121 to 0.0129 all the way to t = 95.
            #
            # **The scaling is the assertion.** ω_B = √(kE₀) = √α here, so the
            # arrest time should go as α^(-1/2), and the accumulated phase
            # ω_B·t at which it happens should not move at all. Measured:
            #
            #   α       0.005   0.01    0.02    0.04
            #   t₀      113.1   73.4    50.5    35.4
            #   t₀√α    8.00    7.34    7.14    7.09
            #
            # A factor of eight in amplitude moves the arrest by 3.2, where
            # α^(-1/2) predicts 2.83. The residue is not something the decay
            # factor in `trapping_phase` removes: evaluated at these arrest
            # times it reads 4.26, 4.79, 5.28 and 5.72, a spread of 34% against
            # the 13% of the undamped column above. The undamped phase is the
            # better invariant here, and it is the one asserted.
            function arrest_time(α; k = 0.3, tmax)
                L = 2*(2π/k)
                Nx = 108
                Δx = L/Nx
                x = collect(range(Δx; step = Δx, length = Nx))
                v = collect(-6:0.1:6)
                t = collect(0.0:0.05:tmax)
                f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + α*cos(k*x)))'
                ε_e = vlasov_poisson(x, v, f₀, t).ε_e
                p = local_extrema(t, ε_e; tmin = 5.0, tmax = tmax - 1, maxima = true)
                ts, rates = Float64[], Float64[]
                for j in 3:2:length(p)
                    i0, i1 = p[j-2], p[j]
                    push!(ts, 0.5*(t[i0] + t[i1]))
                    push!(rates, -log(ε_e[i1]/ε_e[i0])/(2*(t[i1] - t[i0])))
                end
                j = findfirst(≤(0.0), rates)
                j === nothing && error("the damping never stopped at α = $α within t ≤ $tmax")
                # linear interpolation between the last positive rate and the first
                # non-positive one, so the answer is not quantised to the spacing
                # of the maxima
                return ts[j-1] + (ts[j] - ts[j-1])*rates[j-1]/(rates[j-1] - rates[j]), rates
            end

            scan = [(α, first(arrest_time(α; tmax = tmax)))
                    for (α, tmax) in ((0.005, 130.0), (0.01, 90.0),
                                      (0.02, 65.0), (0.04, 45.0))]
            for (α, t₀) in scan
                println("  α = ", rpad(α, 6), " damping stops at t = ", rpad(round(t₀; digits = 1), 6),
                        "  ω_B·t₀ = ", round(t₀*sqrt(α); digits = 2))
                # The bounce phase at the arrest is the invariant statement.
                @test 6.5 < t₀*sqrt(α) < 8.5
            end

            # And the calibration of the guard the linear cases use, in the
            # guard's own units: `trapping_phase` at every arrest time is well
            # above the 2 the linear runs are held under. Measured minimum 4.26.
            γ₃ = landau_rate(0.3)[1]
            @test minimum(trapping_phase(α, γ₃, t₀) for (α, t₀) in scan) > 2.0

            # And the power law itself, fitted rather than eyeballed: log t₀
            # against log α has slope -0.56 over this range, against the -0.5 a
            # constant ω_B would give.
            X = [log(α) for (α, _) in scan]
            Y = [log(t₀) for (_, t₀) in scan]
            n = length(X)
            slope = (n*sum(X.*Y) - sum(X)*sum(Y))/(n*sum(X.^2) - sum(X)^2)
            println("  fitted d(log t₀)/d(log α) = ", round(slope; digits = 3),
                    "  (α^(-1/2) would give -0.5)")
            @test -0.7 < slope < -0.45

            # The linear runs are on the other side of this: at α = 1e-3 and the
            # k = 0.3 window the bounce phase is 1.70, where in the same units
            # the departure sets in at about 2.8 and the arrest at 4.3 to 5.7.
            @test trapping_phase(1e-3, γ₃, 90.0) < 2.0
            @test trapping_phase(1e-2, γ₃, 50.0) > 3.0
        end

        @testset "Each mode recurs at its own 2π/(kΔv)" begin
            # Recurrence is a property of Δv, and `test_free_streaming.jl`
            # measures it exactly -- with no field. Here it is measured in the
            # self-consistent run, mode by mode, because the total `ε_e` cannot
            # tell two modes apart and that has misled this repository before:
            # the Landau notebook shows ε_e rising again at t ≈ 62 and explains
            # it as the seeded mode returning at "π/(kΔv)". It is not. The
            # seeded k = 0.5 mode returns at 2π/(kΔv) = 125.7, twice as late;
            # what arrives at 62 is the *second harmonic*, which the run
            # generates nonlinearly and which recurs at 2π/(2kΔv) = 62.8.
            #
            # Measured on the notebook's own grid, to t = 140:
            #
            #   mode     |E| at t=0   least      peak       at t     T_R
            #   k = 0.5  1.99e-2      8.7e-9     1.03e-2    128.6    125.7
            #   k = 1.0  2.63e-17     --         7.60e-5    64.3     62.8
            #
            # ("least" is the smallest |E| between the decay and the return.)
            # The harmonic starts at round-off -- it is not seeded -- and around
            # its recurrence it is 21 times the seeded mode, which around the
            # seeded mode's own recurrence is 121 times it. In `ε_e` both are
            # the same bump.
            x = collect(π/8:π/8:8π)
            v = collect(-4:0.1:4)
            t = collect(0.0:0.1:140.0)
            k = 0.5
            f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(k*x)))'
            r = vlasov_poisson(x, v, f₀, t; modes = (k, 2k))
            A₁, A₂ = abs.(r.E_modes[:, 1]), abs.(r.E_modes[:, 2])

            for (label, A, T_R) in (("k = 0.5", A₁, 2π/(k*0.1)),
                                    ("k = 1.0", A₂, 2π/(2k*0.1)))
                lo = findfirst(≥(0.6*T_R), t)
                hi = something(findfirst(≥(1.25*T_R), t), length(t))
                i = argmax(view(A, lo:hi)) + lo - 1
                println("  ", label, ": |E| peaks at t = ", round(t[i]; digits = 1),
                        " against T_R = ", round(T_R; digits = 1),
                        "  (", round(A[i]; sigdigits = 3), " from ",
                        round(A[1]; sigdigits = 3), ")")
                # Within one plasma period of the grid's own recurrence time.
                # The lag is real and physical -- the returning ballistic term
                # has to drive the field back up, which takes a fraction of a
                # period -- so it is bounded rather than asserted away.
                @test abs(t[i] - T_R) < 2π/1.4
            end

            # The seeded mode comes back with half the amplitude it started
            # with, which is the scheme's dissipation over 1257 steps rather
            # than anything about the recurrence. Measured 0.520.
            @test A₁[argmax(A₁[findfirst(≥(100.0), t):end]) +
                      findfirst(≥(100.0), t) - 1]/A₁[1] > 0.4

            # And the statement the notebook got wrong: at the time of the first
            # rise in ε_e, it is the harmonic and not the seeded mode.
            #
            # Compared over an envelope rather than at an instant, deliberately.
            # Each mode oscillates at its own frequency and passes through deep
            # nulls, so a point sample of the ratio swings between 20 and 520
            # across three time units without anything physical changing. The
            # peak over the window is the quantity the recurrence is about.
            # Measured over t ∈ [58, 72]: 3.68e-6 against 7.60e-5, a factor 21.
            lo, hi = findfirst(≥(58.0), t), findfirst(≥(72.0), t)
            peak₁, peak₂ = maximum(A₁[lo:hi]), maximum(A₂[lo:hi])
            println("  over t ∈ [58, 72]: max|E_k| = ", round(peak₁; sigdigits = 3),
                    ", max|E_2k| = ", round(peak₂; sigdigits = 3),
                    ", ratio ", round(peak₂/peak₁; digits = 1))
            @test peak₂ > 10*peak₁

            # The other way round at the seeded mode's own recurrence, which is
            # what makes the pair a statement rather than an observation.
            # Measured over t ∈ [120, 135]: 1.03e-2 against 8.51e-5, a factor 121.
            lo₂, hi₂ = findfirst(≥(120.0), t), findfirst(≥(135.0), t)
            @test maximum(A₁[lo₂:hi₂]) > 10*maximum(A₂[lo₂:hi₂])
        end

        @testset "The plasma echo, with the plasma's own field" begin
            # `test_echo.jl` measures the echo with the field off, against a
            # closed form exact in both amplitudes. This is the same experiment
            # in the self-consistent run, against `echo_second_order`: the seed's
            # filament screened at k₁, the kick screened at k₂, and the echo's
            # own density polarising the plasma at k₃ -- three linear responses
            # of a Maxwellian, composed, with nothing taken from the solver. It
            # is the one comparison in this suite with a kinetic theory beyond
            # linear order, and the benchmark Galeotti, Califano and Pegoraro
            # (2006) proposed for Vlasov codes.
            #
            # The field does not correct the echo, it remakes it. Measured at
            # Nx = 128, Δv = 0.05, Δt = 0.02, α = ε = 0.01, τ = 10:
            #
            #                             peak |A₃|   at t
            #   closed form, field off    5.018e-4    30.19
            #   second-order theory       2.385e-4    29.05
            #   run                       2.374e-4    29.05
            #
            # Half the size, a time unit early, and then it rings: the k₃ mode
            # the echo builds is a Langmuir wave, which goes on oscillating and
            # Landau-damping after its ballistic source has phase-mixed away.
            # The comparison is pointwise over all of τ < t ≤ 40, so the ringing
            # is held to the same standard as the peak.
            fine = self_consistent_echo()
            coarse = self_consistent_echo(Nx = 64)

            function against_theory(r)
                after = findall(>(10.0), r.t)
                t = range(r.t[after[1]], r.t[after[end]]; length = length(after))
                A = r.modes[after, 3]
                setup = (α = 0.01, ε = 0.01, τ = 10.0, k₁ = r.k[1], k₂ = r.k[2])
                theory = echo_second_order(t; setup...).echo
                closed = [echo_closed_form(u; setup...) for u in t]
                i, j, c = argmax(abs.(A)), argmax(abs.(theory)), argmax(abs.(closed))
                return (err = maximum(abs, A .- theory)/abs(theory[j]),
                        size = abs(A[i])/abs(theory[j]), t_run = t[i], t_theory = t[j],
                        vs_closed = abs(A[i])/abs(closed[c]), t_closed = t[c],
                        off = maximum(abs, A .- closed)/abs(closed[c]))
            end
            fr, cr = against_theory(fine), against_theory(coarse)
            for (label, r) in (("Nx = 128", fr), ("Nx =  64", cr))
                println("  ", label, ": pointwise ", round(r.err; sigdigits = 3),
                        " of the peak; peak ", round(r.size; digits = 4),
                        " of the theory's, at t = ", round(r.t_run; digits = 2),
                        " (theory ", round(r.t_theory; digits = 2), ")")
            end
            println("  field off: peak ", round(1/fr.vs_closed; digits = 2), "x the run's, at t = ",
                    round(fr.t_closed; digits = 2), "; pointwise ", round(fr.off; sigdigits = 3))

            # 4.8e-3 of the peak pointwise; the peak 0.47% under the theory's,
            # on the same sample.
            @test fr.err < 1e-2
            @test abs(fr.size - 1) < 1e-2
            @test abs(fr.t_run - fr.t_theory) ≤ 2*0.02

            # And it is the field that does it. The closed form has everything
            # else -- the same seed, kick, filament and phase mixing -- and puts
            # the peak at 2.1 times the run's and 1.14 later, missing the run
            # pointwise by 1.01 of its own peak. A run that lost the field would
            # sit on the closed form and fail the comparison above; these hold
            # the two curves apart, so that no run can pass both.
            @test fr.vs_closed < 0.6
            @test fr.t_closed - fr.t_run > 0.5
            @test fr.off > 0.5

            # The residual is the run's truncation, not the theory's: 1.61e-2 at
            # Nx = 64 against 4.81e-3 at 128, a factor 3.4. At Nx = 64, halving Δv
            # or Δt instead leaves it at 1.53e-2 and 1.65e-2, cutting α tenfold
            # moves it by 3e-5, and halving ε by 8e-4 -- the curvature of J₁.
            @test cr.err > 2.5*fr.err
        end

        @testset "A nonlinear equilibrium stays put, and gives on the separatrix" begin
            # Every other run here starts away from equilibrium and is judged on
            # how it moves. This one starts *on* a nonlinear equilibrium -- a
            # function of the particle energy, held by the ion background built
            # for it in `bgk_equilibrium` -- and is judged on how little it moves.
            # It is also the first run in the suite whose physics is carried by
            # trapped particles: at ψ = 0.5 two thirds of them are, and t = 50 is
            # 2.8 periods of the deepest.
            field(r) = maximum(abs, r.E .- r.E₀)/abs(r.E₀)
            drift(r) = maximum(abs, r.f .- r.f₀)/maximum(r.f₀)
            function off_separatrix(r)
                I = argmax(abs.(r.f .- r.f₀))
                return abs(abs(r.v[I[1]]) - r.v_sep[I[2]])/(r.v[2] - r.v[1])
            end
            fine = (Nx = 128, Δv = 0.05, Δt = 0.025)
            smooth, smooth_fine = bgk_equilibrium(), bgk_equilibrium(; fine...)
            kinked = bgk_equilibrium(T_trapped = 2.0)
            kinked_fine = bgk_equilibrium(; T_trapped = 2.0, fine...)
            for (label, r) in (("Maxwell-Boltzmann, 64 x 121  ", smooth),
                               ("Maxwell-Boltzmann, 128 x 241 ", smooth_fine),
                               ("trapped at T = 2, 64 x 121   ", kinked),
                               ("trapped at T = 2, 128 x 241  ", kinked_fine))
                println("  ", label, " f ", rpad(round(drift(r); sigdigits = 3), 10),
                        " field ", rpad(round(field(r); sigdigits = 3), 10),
                        " worst cell ", round(off_separatrix(r); digits = 1),
                        " cells off the separatrix")
            end

            @testset "the smooth one holds, and converges at the scheme's order" begin
                # Maxwell–Boltzmann, analytic across the separatrix. Through
                # t = 50: f within 1.52e-3 of its peak, the field within 3.85e-3.
                # Neither is an oscillation about the equilibrium; both grow
                # steadily -- the field's departure is 6.5e-3 by t = 100 -- and
                # they are the scheme's dissipation, since L² falls by 1.13e-3
                # and the entropy rises by 1.27e-4 where the exact flow keeps both.
                @test drift(smooth) < 5e-3
                @test field(smooth) < 1e-2
                @test smooth.l2[end] < smooth.l2[1]
                @test smooth.entropy[end] > smooth.entropy[1]

                # Halving Δx, Δv and Δt together: 7.6 times less in f and 6.1 in
                # the field, PFC's third order. Halving Δt alone moves neither
                # (3.85e-3 → 3.89e-3), so the splitting is not what limits it.
                @test drift(smooth)/drift(smooth_fine) > 5
                @test field(smooth)/field(smooth_fine) > 4
            end

            @testset "a kink on the separatrix is where it gives" begin
                # Trapped particles at twice the temperature of the passing ones.
                # F is continuous and still decreasing -- so the equilibrium is
                # still stable -- but its slope jumps where trapped meets passing,
                # and that is where the error goes: the worst cell is on the
                # separatrix at both resolutions, where the smooth equilibrium's
                # is 13 and 24 cells away from it.
                @test off_separatrix(kinked) ≤ 1
                @test off_separatrix(kinked_fine) ≤ 1
                @test off_separatrix(smooth) > 5

                # And it converges at first order, not third: 1.90e-2 → 1.09e-2,
                # a factor 1.75, twelve times the smooth case's error to begin
                # with. The field hardly notices, 3.69e-3 → 1.43e-3: it is an
                # integral over the kink, not a sample of it.
                @test 1.3 < drift(kinked)/drift(kinked_fine) < 3
                @test drift(kinked) > 5*drift(smooth)
                @test field(kinked) < 1e-2
            end

            @testset "a trapped population above f = 1 holds as well" begin
                # Trapped particles at half the temperature peak at 1.79, above
                # the bound of 1 the harness used to give PFC. Nothing checked
                # it, and the limiter wrecked the run: 43.6% of the peak by
                # t = 50. With the bound taken from f₀, as it is now: 0.99%, on
                # the separatrix again.
                cold = bgk_equilibrium(T_trapped = 0.5)
                @test maximum(cold.f₀) > 1.5
                @test drift(cold) < 3e-2
                @test off_separatrix(cold) ≤ 1
            end

            @testset "and it is this equilibrium, not any" begin
                # Ions built with the Poisson sign `docs/normalization.md` used
                # to give, ∂E/∂x = nᵢ − nₑ: the code's field is then the
                # equilibrium's reversed from the first sample (E/E₀ = −1.000),
                # and f departs by 16% of its peak by t = 50.
                reversed = bgk_equilibrium(ion_sign = -1)
                @test field(reversed) > 1.5
                @test drift(reversed) > 0.1

                # A potential 10% off the one the ions hold puts the field 49% and
                # 59% off at once -- the field is the small difference nₑ − nᵢ --
                # and f departs by 3.7% and 4.0%, 24 and 27 times the
                # equilibrium's own drift.
                for ψ in (0.45, 0.55)
                    r = bgk_equilibrium(; ψ, ψ_ions = 0.5)
                    @test field(r) > 0.3
                    @test drift(r) > 10*drift(smooth)
                end
            end

            @testset "the ions it is handed are the ions it keeps" begin
                # Handed `nᵢ`, the driver takes `f` as matched to it. Its
                # rescaling to the ions' charge is a trapezoid, exact only for
                # proportional profiles, and on this pair it is 1 − 1.9e-3:
                # applied, it doubles the drift above, 3.85e-3 → 8.07e-3 in the
                # field and 1.52e-3 → 3.17e-3 in f, and both stay inside their
                # tolerances. So it is pinned here, one step in, on the flux-form
                # mass: f₀'s to 2.2e-16 by default, 1.9e-3 short when forced.
                one_step(; kw...) = vlasov_poisson(smooth.x, smooth.v, smooth.f₀, [0.0, 0.05];
                                                   smooth.nᵢ, invariants = true, kw...).mass[1]
                cells = (smooth.v[2] - smooth.v[1])*(smooth.x[2] - smooth.x[1])
                Σf₀ = sum(smooth.f₀)*cells
                @test one_step() ≈ Σf₀ rtol = 1e-12
                @test one_step(renormalize = true) < (1 - 1e-3)*Σf₀
            end
        end

        @testset "The Vlasov-Poisson flow is reversible, and what breaks it" begin
            # `test_strang_splitting.jl` measures reversibility on a rigid
            # rotation with the field switched off. This is the same statement
            # for the self-consistent system: Vlasov--Poisson is invariant under
            # `(t, v) → (-t, -v)` with `E` unchanged, so running forward,
            # flipping the velocity axis, running forward again and flipping
            # back must return the initial state. Strang splitting is symmetric
            # and preserves that exactly; what does not is the scheme's own
            # dissipation, which has no time reverse.
            #
            # The flip is `f[end:-1:1, :]` on a velocity grid symmetric about
            # zero -- asserted below, since on any other grid it would be off by
            # a fraction of a cell and the test would measure that instead.
            k = 0.5
            L = 2*(2π/k)
            function round_trip(Nx, Δv, Δt, T, α; scheme = nothing)
                Δx = L/Nx
                x = collect(range(Δx; step = Δx, length = Nx))
                v = collect(-6:Δv:6)
                @assert v[1] == -v[end] && iseven(length(v) - 1)
                t = collect(0.0:Δt:T)
                f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + α*cos(k*x)))'
                flip(f) = f[end:-1:1, :]
                fwd = vlasov_poisson(x, v, f₀, t; scheme_x = scheme, scheme_v = scheme,
                                     invariants = true)
                back = flip(vlasov_poisson(x, v, flip(fwd.f), t;
                                           scheme_x = scheme, scheme_v = scheme).f)
                # the driver renormalises what it is handed, so the comparison is
                # against `f₀` at the normalisation the round trip came back with
                ref = f₀ .* (sum(back)/sum(f₀))
                return (; err = maximum(abs, back .- ref)/maximum(f₀),
                          fmin = minimum(fwd.fmin[1:end-1]),
                          l2 = (fwd.l2[end-1] - fwd.l2[1])/fwd.l2[1])
            end

            @testset "at small amplitude it converges away" begin
                # Nothing is irreversible about the equations, so the round-trip
                # error is the scheme's dissipation and has to vanish with the
                # grid. Measured at α = 0.05 over T = 20: 2.82e-3, 6.48e-4 and
                # 1.58e-4 at Nx = 64, 128 and 256 -- ratios 4.4 and 4.1, second
                # order. The last level costs 22 s and is left out; the two below
                # cost two seconds.
                #
                # `PFC` is third order, and with its bound at 1, which never
                # engaged, this read 1.97e-3, 3.44e-4 and 4.85e-5 -- ratios 5.7
                # and 7.1. The bound is now the initial maximum, as Liouville
                # has it, and the limiter clips the peak cell there; the error is
                # a maximum over f, so it sees that cell, and one order is what
                # the clipping costs. The threshold is set for second order.
                coarse = round_trip(64, 0.1, 0.05, 20.0, 0.05)
                fine = round_trip(128, 0.05, 0.025, 20.0, 0.05)
                println("  α = 0.05: round trip ", round(coarse.err; sigdigits = 3),
                        " -> ", round(fine.err; sigdigits = 3),
                        "  (x", round(coarse.err/fine.err; digits = 1),
                        "; second order gives 4, third 8)")
                @test fine.err < coarse.err/3.5
                @test fine.err < 1e-3
            end

            @testset "at large amplitude the grid loses it, and refining does not help" begin
                # The same measurement at α = 0.5, where the flow folds the
                # distribution into filaments finer than Δv within a few plasma
                # periods. Past that the information needed to run the film
                # backwards is not on the grid any more, and the round trip
                # returns 16% of the peak whatever the resolution: measured
                # 1.64e-1 at Nx = 64 against 1.56e-1 at 128, a factor of 1.05
                # where the linear case gains 4.4.
                #
                # This is the honest counterweight to the testset above. The
                # scheme is third order and the splitting is symmetric, and
                # neither buys reversibility in a run that has made structure
                # below the mesh -- which is what every nonlinear run in this
                # suite is doing by the time it is interesting.
                coarse = round_trip(64, 0.1, 0.05, 20.0, 0.5)
                fine = round_trip(128, 0.05, 0.025, 20.0, 0.5)
                println("  α = 0.5:  round trip ", round(coarse.err; sigdigits = 3),
                        " -> ", round(fine.err; sigdigits = 3),
                        "  (x", round(coarse.err/fine.err; digits = 2), ")")
                @test coarse.err > 0.1
                @test coarse.err/fine.err < 1.5
            end

            @testset "and the schemes that keep it are the ones that go negative" begin
                # Ranked by round-trip error at α = 0.5, with what each scheme
                # costs to get there:
                #
                #   scheme                 round trip   min f       ΔL²/L²
                #   LaxWendroff            8.50e-2      -9.44e-2    -0.005
                #   SemiLagrangian cubic   1.06e-1      -5.82e-2    -0.004
                #   PFC                    1.64e-1      +5.21e-10   -0.048
                #   Upwind                 3.32e-1      +3.10e-09   -0.210
                #
                # The ordering is the one `verification/scheme-comparison.jl`
                # reports on the damping rate, arrived at through a completely
                # different quantity: the two schemes that best preserve the
                # flow are the two that are not monotone, and they pay for it by
                # driving `f` to -16% of its peak. Reversibility and positivity
                # are the same trade-off seen from two sides, and `PFC` is the
                # default because a negative distribution has no entropy.
                results = Dict{String,Any}()
                for (name, scheme) in (("LaxWendroff", LaxWendroff()),
                                       ("SemiLagrangian cubic", SemiLagrangian(CubicSpline())),
                                       ("PFC", f -> PFC(fmin = 0.0, fmax = maximum(f))),
                                       ("Upwind", Upwind()))
                    r = round_trip(64, 0.1, 0.05, 20.0, 0.5; scheme = scheme)
                    results[name] = r
                    println("  ", rpad(name, 22), "round trip ", rpad(round(r.err; sigdigits = 3), 9),
                            " min f = ", rpad(round(r.fmin; sigdigits = 3), 11),
                            " ΔL²/L² = ", round(r.l2; digits = 3))
                end
                @test results["LaxWendroff"].err < results["PFC"].err
                @test results["SemiLagrangian cubic"].err < results["PFC"].err
                @test results["Upwind"].err > 2*results["PFC"].err
                # and the price, which is why the default is the slower one
                @test results["LaxWendroff"].fmin < -0.05
                @test results["SemiLagrangian cubic"].fmin < -0.05
                @test results["PFC"].fmin ≥ 0.0
                # upwind keeps positivity and loses the physics instead: a fifth
                # of the L² norm, four times what PFC dissipates
                @test results["Upwind"].l2 < 4*results["PFC"].l2
            end
        end

        @testset "Strong Landau damping: the damping stops and reverses" begin
            # The `α = 0.5` case, which the notebook has run since 2021 and
            # compared "by eye" against Fig. 6(a) of Filbet, Sonnendrücker and
            # Bertrand (2001). It is the standard nonlinear benchmark -- the
            # field damps, the resonant particles trap, and the field grows
            # again -- and it is quoted in the literature by two numbers, so
            # there was no reason for it to stay a picture.
            #
            # Reported values scatter with the window and the resolution, which
            # this testset is in a position to show rather than gloss:
            # Cheng and Knorr (1976) give γ₁ = -0.281 and γ₂ = 0.084, and later
            # work quotes -0.292 with 0.0815 and -0.2918 with 0.08584.
            #
            # Measured here, on 128 x 241 with Δt = 0.025:
            #
            #   γ₁ over [0.5, 12]   0.2863     (4 maxima)
            #   γ₂ over [20, 40]    0.0789     (8 maxima)
            #
            # **Both windows are conventions, and the sensitivity is the reason
            # to say so.** γ₁ reads 0.3786 over a window holding three maxima
            # and 0.2281 over one holding five: the envelope is not an
            # exponential, it steepens and then flattens into the trapping
            # plateau, so a "damping rate" is a straight line through a curve
            # and the answer depends on how much of the curve is in the window.
            # The window here is the one that spans the decay proper -- four
            # maxima, ending before the plateau at t ≈ 13 -- and it puts γ₁ in
            # the literature's range rather than beside it. γ₂ is not in the
            # range the citations above give: 0.0789 sits 3.2% under the lowest
            # of them, 0.0815, and it is refinement rather than the window that
            # closes the gap (below).
            #
            # γ₂ is better behaved (0.0751 to 0.0789 over windows holding eight
            # or nine maxima) but saturates after t ≈ 41, where the field stops
            # growing: [20, 44] reads 0.0657 because it averages the turnover
            # in.
            k = 0.5
            L = 2*(2π/k)
            function strong_case(Nx, Δv, Δt; v = collect(-6:Δv:6))
                Δx = L/Nx
                x = collect(range(Δx; step = Δx, length = Nx))
                t = collect(0.0:Δt:45.0)
                f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.5*cos(k*x)))'
                r = vlasov_poisson(x, v, f₀, t; invariants = true)
                γ₁, n₁ = damping_rate(t, r.ε_e; tmin = 0.5, tmax = 12.0)
                γ₂, n₂ = damping_rate(t, r.ε_e; tmin = 20.0, tmax = 40.0)
                return (; γ₁, n₁, γ₂ = -γ₂, n₂, r)
            end

            fine = strong_case(128, 0.05, 0.025)
            println("  128 x 241: γ₁ = ", round(fine.γ₁; digits = 4), " (", fine.n₁,
                    " maxima), γ₂ = ", round(fine.γ₂; digits = 4), " (", fine.n₂, " maxima)")
            @test 0.27 < fine.γ₁ < 0.30       # the literature's -0.281 to -0.292
            @test 0.070 < fine.γ₂ < 0.090     # cited 0.0815 to 0.08584; 0.0789 here

            # Refinement moves γ₂ toward the published value rather than away
            # from it, which is the statement that the agreement is not a
            # coincidence of this grid: 0.0716 at half the resolution, 0.0789
            # here, 0.0814 at twice it (that last run costs 22 s and is not
            # repeated in CI). γ₁ is converged already -- 0.2794 against 0.2863
            # -- because the first decade of the decay is resolved on both.
            coarse = strong_case(64, 0.1, 0.05)
            println("  64 x 121:  γ₁ = ", round(coarse.γ₁; digits = 4),
                    ", γ₂ = ", round(coarse.γ₂; digits = 4),
                    "   (γ₂ at twice the fine resolution: 0.0814)")
            @test coarse.γ₂ < fine.γ₂
            @test abs(fine.γ₂ - 0.0815) < abs(coarse.γ₂ - 0.0815)

            # At this amplitude the distribution comes close to zero, which is
            # what `PFC` is for and what the comparison study measures the other
            # schemes failing. Measured minimum: 1.9e-9, and mass to round-off.
            println("  min f = ", minimum(fine.r.fmin[1:end-1]),
                    "   mass drift = ",
                    maximum(abs, fine.r.mass[1:end-1] .- fine.r.mass[1])/fine.r.mass[1])
            @test minimum(fine.r.fmin[1:end-1]) ≥ 0.0
            @test maximum(abs, fine.r.mass[1:end-1] .- fine.r.mass[1])/fine.r.mass[1] < 1e-13

            # And the non-uniform velocity grid -- the one the notebook runs,
            # and the only path this suite has toward an adaptive mesh -- gives
            # the same physics. Until now it was asserted only through an energy
            # drift over a *linear* run, where the distribution never approaches
            # the sharp gradients its limiter exists for. Measured against the
            # uniform grid at the same Δt: γ₁ 0.2793 against 0.2794 (0.05%) and
            # γ₂ 0.0721 against 0.0716 (0.6%).
            #
            # This grid is the one run in the suite whose velocity sweep starts
            # past its Courant limit: the field's amplitude is 1.0017 on the
            # first step and Δt = 0.05 is the width of the narrow cells, so it
            # asks for 1.0017 of them. `advect!` refuses that, `line_advector`
            # splits those four calls in two, and γ₁ and γ₂ move in the sixth
            # digit (0.279257 → 0.279258, 0.072051 → 0.072050).
            stretched = strong_case(64, 0.1, 0.05;
                v = vcat(collect(-6:0.1:-1.1), collect(-1:0.05:1), collect(1.1:0.1:6)))
            println("  non-uniform Δv: γ₁ = ", round(stretched.γ₁; digits = 4),
                    ", γ₂ = ", round(stretched.γ₂; digits = 4),
                    "   against the uniform grid at the same Δt")
            @test isapprox(stretched.γ₁, coarse.γ₁; rtol = 0.02)
            @test isapprox(stretched.γ₂, coarse.γ₂; rtol = 0.02)
        end

        @testset "A drifting plasma damps the same way, Doppler-shifted" begin
            # Every Vlasov--Poisson case in this suite starts from a
            # distribution symmetric in `v`, so `u = 0` throughout and the
            # drifting half of the solver is never exercised. That is exactly
            # the blind spot `test_damping_1v.jl` found in `BGK`, where the
            # mean-velocity computation had never run on data with a mean
            # velocity; here it would hide an error in the `v` sweep that
            # cancels between `+v` and `-v`, or a resonance found by symmetry
            # rather than by physics.
            #
            # The statement is Galilean invariance. Boosting by `u` carries
            #
            #     f(x, v, t) → f(x - ut, v - u, t),   E(x, t) → E(x - ut, t)
            #
            # so the field is the same solution translated: `|E_k|` is
            # unchanged, the damping rate is unchanged, and the only difference
            # is a phase `exp(-ikut)` on the mode. Nothing in the *discretisation*
            # is Galilean invariant -- the grid does not move, and the boosted
            # Maxwellian sits on it asymmetrically -- so the agreement below is
            # a measurement rather than an identity.
            #
            # Measured at the 11 maxima of `|E_k|` in the fitting window, with
            # `u = 0.5` against `u = 0`:
            #
            #   amplitude ratio     within 1.7e-3
            #   phase               within 1.8e-3 rad
            #   fitted γ            0.136% apart
            #   fitted ω            identical to the estimator's resolution
            #
            # **The Doppler factor is removed at the time the field was sampled,
            # `t[k] + Δt/2`**, not at `t[k]`: `E_modes` is recorded from the field
            # solved mid-step (see `vlasov_poisson`). Removed at `t[k]`, the same
            # comparison reads 8.0e-3 rad, and this comment once quoted that as
            # the grid's non-invariance -- 78% of it was `k·u·Δt/2 = 6.25e-3`, a
            # constant offset from the timestamp. The tolerance is set against
            # the corrected figure, so an offset of that size now fails.
            #
            # Compared at the maxima, and deliberately: both `|E_k|` and `ε_e`
            # pass through deep nulls, where a relative difference of anything
            # is meaningless. The first attempt at this test compared them
            # pointwise and reported 900% -- entirely from two nulls landing a
            # time step apart.
            k, α, u = 0.5, 1e-3, 0.5
            L = 2*(2π/k)
            Nx = 64
            Δx = L/Nx
            x = collect(range(Δx; step = Δx, length = Nx))
            v = collect(-6:0.1:6)       # wide enough for the boosted resonance at 3.33
            t = collect(0.0:0.05:40.0)
            drifting(u) = vlasov_poisson(x, v,
                1/sqrt(2π)*(@. exp(-0.5*(v - u)^2)) * (@. (1.0 + α*cos(k*x)))',
                t; modes = (k,))
            rest, boosted = drifting(0.0), drifting(u)

            # The Doppler factor is taken out here, at each sample's own time;
            # the assertion is that what remains is the same complex history.
            Δt = t[2] - t[1]
            A₀ = rest.E_modes[:, 1]
            A_u = boosted.E_modes[:, 1] .* cis.(k*u .* (t .+ Δt/2))
            peaks = local_extrema(t, abs.(A₀); tmin = 5.0, tmax = 30.0, maxima = true)
            amp = maximum(abs(abs(A_u[i])/abs(A₀[i]) - 1) for i in peaks)
            phase = maximum(abs(angle(A_u[i]/A₀[i])) for i in peaks)
            println("  boosted by u = ", u, ", over ", length(peaks), " maxima: amplitude ",
                    round(amp; sigdigits = 3), ", phase ", round(phase; sigdigits = 3), " rad")
            @test amp < 5e-3
            @test phase < 5e-3

            # And the rates, which is the same statement read through the
            # estimators the rest of this file uses.
            γ₀, _ = damping_rate(t, rest.ε_e; tmin = 6.0, tmax = 30.0)
            γ_u, _ = damping_rate(t, boosted.ε_e; tmin = 6.0, tmax = 30.0)
            ω₀, _ = oscillation_frequency(t, rest.ε_e; tmin = 6.0, tmax = 30.0)
            ω_u, _ = oscillation_frequency(t, boosted.ε_e; tmin = 6.0, tmax = 30.0)
            println("  γ = ", round(γ₀; digits = 5), " at rest, ", round(γ_u; digits = 5),
                    " boosted (", round(100*abs(γ_u - γ₀)/γ₀; digits = 3), "%)",
                    "   ω = ", round(ω₀; digits = 5), " and ", round(ω_u; digits = 5))
            @test isapprox(γ_u, γ₀; rtol = 5e-3)
            @test isapprox(ω_u, ω₀; rtol = 1e-3)

            # Teeth: the Doppler shift being taken out is a real shift, not a
            # formality. Without the correction the mode's phase differs by up
            # to 2.88 radians over the same window.
            raw = maximum(abs(angle(boosted.E_modes[i, 1]/A₀[i])) for i in peaks)
            println("  without the exp(ikut) correction the phase differs by up to ",
                    round(raw; digits = 2), " rad")
            @test raw > 1.0
        end

        @testset "Landau damping converges under refinement" begin
            # Agreement at one resolution inside a 3% band can be luck: two
            # errors of opposite sign meeting in the middle is exactly how a
            # plausible-but-wrong solver survives a tolerance. What cannot be
            # luck is the error *shrinking* when the grid is refined.
            #
            # **This one keeps α = 1e-2 where the cases above moved to 1e-3**,
            # and the bounce phase says why it may: at k = 0.5 the damping is
            # fast enough that `ω_B` falls with the field before a bounce
            # completes, so `trapping_phase(1e-2, γ, 30) = 0.65` against the 3.7
            # the k = 0.3 case reached. The larger amplitude keeps the mode
            # above the recurrence floor at the coarsest grid, where α = 1e-3
            # would put it into the noise at Nx = 32.
            #
            # Δx, Δv and Δt are halved together, so the Courant number stays at
            # 0.81 and only the discretisation moves. Measured at k = 0.5:
            #
            #   Nx    Δv      Δt     γ         error    ΔL²/L² over the run
            #   32    0.2     0.16   0.16318   6.41%    -1.11e-4
            #   64    0.1     0.08   0.15536   1.31%    -1.01e-5
            #   128   0.05    0.04   0.15430   0.61%    -1.27e-6
            #   256   0.025   0.02   0.15407   0.47%    -1.63e-7
            #
            # The last level costs 14 s on its own and is left out; the three
            # below cost about two seconds together.
            #
            # The L² column is why the γ column behaves as it does, and is worth
            # asserting alongside it. `PFC` is third order, so halving the grid
            # should cut its dissipation by eight -- measured 11.0, 7.9 and 7.8.
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
                x = collect(range(Δx; step = Δx, length = Nx))
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
            # reports 2.4e-4 of mass drift that belongs to the quadrature rather
            # than the solver. See the note on `vlasov_poisson`.
            #
            # Measured over the k = 0.5 case, 875 steps:
            #
            #   mass       2.8e-16 relative               -- round-off
            #   momentum   5.3e-16 absolute, on mass 25.1 -- round-off
            #   L²         -1.0e-5, monotone decreasing   -- numerical dissipation
            #   entropy    +7.5e-6, monotone increasing   -- the same thing
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
            Nx = 64
            Δx = L/Nx
            x = collect(range(Δx; step = Δx, length = Nx))
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
            # difference is not cosmetic: at the `vt = 0.3` of the three
            # growth-rate runs the cold form is off by up to 3.14% where the
            # warm root is off by 1.95%, and at the `vt = 0.6` run by 7.44%
            # against 2.01% -- the cold error grows with the temperature, the
            # warm one does not. The one place it changes a *conclusion* rather
            # than a number is `a = 0.8`; see below.

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
                #   0.4       0.30173      0.30362   -0.62%    0.30819   -2.10%
                #   0.6       0.34228      0.34909   -1.95%    0.35339   -3.14%
                #   0.8       0.31229      0.31201   +0.09%    0.31134   +0.31%
                #
                # The a = 0.4 row read 0.30245, -0.39%, while `two_stream` built
                # that grid a point short: 95 cells, whose fundamental is
                # a = 0.4042, and against the root there it was -0.95%.
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
                # the warm root the same measurement is +0.09%, and nothing
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
                # follow it. Measured γ = 0.32710 at vt = 0.6 against 0.34228 at
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
                # a = 1.2 gives 1.00 (4.84e-5 decaying to 2.54e-6) and a = 1.6
                # gives 1.00 (2.02e-5 to 1.2e-9), against 2.0e8 at a = 0.6 over
                # its own t ≤ 24.
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
                # Measured γ = 0.09516 against 0.09823, which is 3.12%, fitted
                # over t ∈ [45.6, 78.2]. The band is the same `[100ε₀, 5.0]`
                # every other case uses; `tmax = 80` is what it takes to reach
                # 5.0 at a tenth of the growth rate, and the velocity sweep stays
                # under its Courant limit throughout, at 0.94 at most.
                #
                # Held to 8% rather than the 3% above, and the reason is in
                # `growth_rate`: the beat between the growing root and the
                # oscillating pair decays as exp(-γt) relative to the mode, so
                # it is worst where γ is smallest, and γ here is a third of the
                # branch maximum. Measured over `hi` ∈ {1, 2, 3, 5} the fit
                # moves from -6.01% to -3.12%; 8% covers that spread with room.
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
