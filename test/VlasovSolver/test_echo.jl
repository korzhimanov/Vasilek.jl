using Vasilek

@isdefined(echo_closed_form) || include(joinpath(@__DIR__, "..", "echo.jl"))

"""
The plasma echo: phase mixing run backwards by a second perturbation.

`test_free_streaming.jl` shows a density mode decaying as `exp(−k²t²/2)` with
nothing dissipated: the perturbation is carried into ever finer structure in `v`,
whose velocity integral cancels. This file shows the information is still there.
A kick `v → v + ε cos k₂x` at `t = τ` shears against the filaments of the seeded
`k₁` mode, and their `k₂ − k₁` beat un-mixes: a density mode nobody seeded
appears around

    t_e = k₂τ/(k₂ − k₁) = 15

and phase-mixes away again, with the amplitude of `echo_closed_form`.

What makes it a test of the solver rather than of free streaming is what has to
survive for the echo to happen. At the kick the seeded mode is down to
`α·exp(−k₁²τ²/2) = 3.7e-7`, where no moment of `f` tells it from zero, and the
kick has to carry filaments of wavelength `2π/(k₁τ) = 1.26` in `v` -- twelve cells
at `Δv = 0.1` -- without smearing them. Numerical diffusion acts on those
filaments as a collision operator would, which is why echoes are the laboratory
measure of weak collisionality (Su & Oberman 1968), and why Galeotti, Califano and
Pegoraro (2006) proposed them as a benchmark for Vlasov codes.

The setup and the closed form are in `test/echo.jl`.
"""

"Pointwise error against the closed form over `t ≥ τ` relative to its peak, and the measured peak."
function echo_error(r; α = 0.1, ε = 0.2, τ = 5.0)
    k₁, k₂, _ = r.k
    exact = [echo_closed_form(s; α, ε, τ, k₁, k₂) for s in r.t]
    after = findall(≥(τ), r.t)
    A = r.modes[:, 3]
    i = after[argmax(abs.(A[after]))]
    return (curve = maximum(abs, A[after] .- exact[after])/maximum(abs, exact),
            peak = A[i], t_peak = r.t[i])
end

@testset "Plasma echo" begin
    # PFC bounded by the distribution it carries, whose maximum is (1 + α)/√(2π);
    # built by the driver from the `f` it actually starts on.
    pfc = f -> PFC(fmin = 0.0, fmax = maximum(f))

    # The production scheme at three velocity resolutions, shared by the first
    # and the last testset. Nx = 128 rather than 64, because at 64 the x-sweep's
    # own truncation (2.8e-3 of the peak) is a floor under the two finest; and
    # Δt = 0.01, because PFC is restricted to Courant ≤ 1 and at Nx = 128 the
    # fastest row would carry 1.22 at Δt = 0.02.
    Δvs = (0.2, 0.1, 0.05)
    elapsed = @elapsed refined = [ballistic_echo(pfc, pfc; Nx = 128, Δt = 0.01, Δv) for Δv in Δvs]
    errors = echo_error.(refined)
    for (Δv, e) in zip(Δvs, errors)
        println("  PFC, Δv = ", rpad(Δv, 5), " pointwise error ", round(e.curve; sigdigits = 3),
                " of the peak, peak at t = ", e.t_peak)
    end
    println("  (", round(elapsed; digits = 2), " s)")

    @testset "the echo comes when and as the closed form says" begin
        k₁, k₂, k₃ = refined[end].k
        dense = range(5.0, 22.0; length = 17001)
        exact = abs.(echo_closed_form.(dense; α = 0.1, ε = 0.2, τ = 5.0, k₁, k₂))
        linear = [0.1*(k₃*0.2*(s - 5.0)/2)*exp(-(k₃*s - k₂*5.0)^2/2) for s in dense]
        t_exact = dense[argmax(exact)]
        fine = errors[end]

        # Measured at Δv = 0.05: 1.25e-3 of the peak, pointwise over t ≥ τ.
        @test fine.curve < 3e-3

        # The peak at t = 15.28, where the closed form's is, and not at t_e = 15,
        # fourteen steps earlier: J₁ is still growing across the pulse. The
        # linear-in-ε form peaks at 15.39.
        @test abs(fine.t_peak - t_exact) ≤ 2*0.01
        @test abs(fine.t_peak - k₂*5.0/k₃) > 10*0.01

        # The Bessel function and not its small-argument limit: J₁(z) → z/2 puts
        # the peak 14.6% high, and the run is within 0.11% of J₁.
        @test abs(abs(fine.peak)/maximum(exact) - 1) < 3e-3
        @test maximum(linear)/maximum(exact) - 1 > 0.1

        # The sign. The echo is +sin k₃x, so its amplitude is negative imaginary --
        # to 4e-15 of its modulus, the setup being mirror-symmetric -- and
        # reversing the kick reverses it. A kick applied as v → v − ε cos k₂x
        # leaves |A₃|, and so the peak's time and size, exactly as they were; the
        # pointwise comparison catches it as an error of twice the peak, and these
        # say which way round it went.
        @test imag(fine.peak) < 0
        @test abs(real(fine.peak)) < 1e-10*abs(fine.peak)
        flipped = echo_error(ballistic_echo(pfc, pfc; ε = -0.2); ε = -0.2)
        @test imag(flipped.peak) > 0
        @test flipped.curve < 3e-2
    end

    @testset "a time the run never reaches is refused" begin
        # A snapshot past `tmax`, or between two steps, used to come back as
        # uninitialised memory; the kick past `tmax` as a BoundsError. Both are
        # caught before anything runs, so these cost nothing.
        @test_throws ErrorException ballistic_echo(pfc, pfc; snapshots = (25.0,))
        @test_throws ErrorException ballistic_echo(pfc, pfc; snapshots = (11.005,))
        @test_throws ErrorException ballistic_echo(pfc, pfc; τ = 30.0)
        @test_throws ErrorException ballistic_echo(pfc, pfc; τ = 5.005)
    end

    @testset "it is made of what no moment could see" begin
        # At the kick the seeded mode's density amplitude is 1.84e-6 -- the closed
        # form says 3.7e-7, and the rest is the x-sweep damping rows unevenly so
        # their phases no longer cancel -- and the echo it turns into is 0.0444,
        # 2.4e4 times larger.
        r = refined[end]
        hidden = abs(r.modes[r.kick, 1])
        println("  seeded mode at the kick ", round(hidden; sigdigits = 3),
                ", echo ", round(abs(errors[end].peak); sigdigits = 3))
        @test hidden < 1e-5
        @test abs(errors[end].peak) > 1e4*hidden
    end

    @testset "what a scheme keeps of a filament is what the echo returns" begin
        # Measured at Nx = 64, Δv = 0.1, pointwise over t ≥ τ against the peak:
        #
        #   SemiLagrangian cubic  1.07e-3
        #   PFC                   1.03e-2
        #   LaxWendroff           4.48e-2
        #   Upwind                0.449    returning 0.553 of the peak
        #
        # The same order as on the free-streaming mode, but not the same
        # distances: upwind's 1.7e-2 there is 45% here. Most of that is its
        # x-sweep, which diffuses each row's spatial modulation at a rate
        # proportional to |v|, so the rows lose phase coherence and not only
        # amplitude -- at the kick it has left 1.8e-2 of the seeded mode in the
        # density, where the exact solution has 3.7e-6. Refining Δv alone takes
        # its loss only from 0.52 to 0.39.
        cases = [("SemiLagrangian cubic", SemiLagrangian(CubicSpline()), 3e-3),
                 ("PFC",                  pfc,                            3e-2),
                 ("LaxWendroff",          LaxWendroff(),                  0.1),
                 ("Upwind",               Upwind(),                       0.6)]
        errs = Dict{String,Float64}()
        for (name, scheme, tol) in cases
            e = echo_error(ballistic_echo(scheme, scheme))
            errs[name] = e.curve
            returned = abs(e.peak)/abs(echo_closed_form(15.28; α = 0.1, ε = 0.2, τ = 5.0,
                                                        k₁ = 1.0, k₂ = 1.5))
            println("  ", rpad(name, 22), "pointwise error ", round(e.curve; sigdigits = 3),
                    ", peak returned ", round(returned; digits = 4))
            @test e.curve < tol
            name == "Upwind" && @test returned < 0.7
        end
        @test errs["SemiLagrangian cubic"] < errs["PFC"] < errs["LaxWendroff"] < errs["Upwind"]
    end

    @testset "and the loss is truncation error, at PFC's order" begin
        # 6.39e-2 → 8.61e-3 → 1.25e-3 as Δv halves: ratios 7.4 and 6.9, orders
        # 2.89 and 2.78 against the scheme's 3. The second ratio is the lower
        # because the x-sweep's floor at Nx = 128 is a quarter of the finest
        # error: refining Δv further, to 0.025 and 0.0125, gives 4.3e-4 and
        # 3.4e-4.
        @test errors[1].curve/errors[2].curve > 5
        @test errors[2].curve/errors[3].curve > 5
    end

    @testset "the self-consistent theory checks out before anything runs against it" begin
        # `echo_second_order` is what `test_verification.jl` holds the
        # self-consistent run to, so its pieces are checked here first, each
        # against something that takes a different route.
        maxwellian = ((density = 1.0, drift = 0.0, vt = 1.0),)
        grid = [(k, v) for k in (0.5, 1.0, 1.5), v in -4:0.25:4]

        # The dielectric function at the resonance is the package's, to
        # 1.2e-16, and its ω-derivative is that function's, to 8.3e-10 against a
        # central difference.
        @test maximum(abs(resonant_dielectric(k, v)[1] - dielectric(k*v, k, maxwellian))
                      for (k, v) in grid) < 1e-12
        @test maximum(abs(resonant_dielectric(k, v)[2] -
                          (dielectric(k*v + 1e-6, k, maxwellian) -
                           dielectric(k*v - 1e-6, k, maxwellian))/2e-6)
                      for (k, v) in grid) < 1e-7

        # The filament's closed form against the general form it came from,
        # (α/2)M − M′·ê₁(k₁v), with the seed's field solved in the time domain by
        # `maxwellian_response` and transformed by quadrature. Measured 4.9e-6,
        # and second order in the step: 1.95e-5 at Δt = 0.01, 7.8e-7 at 0.002.
        # This is what ties the Volterra kernel -- sign, normalisation and all --
        # to the Z-function side.
        Δt = 0.005
        s = 0:Δt:30.0
        e = maxwellian_response([exp(-u^2/2)/2 + 0im for u in s], s, 1.0) ./ im
        worst = 0.0
        for v in (-2.0, -0.7, 0.0, 0.5, 1.3, 2.5)
            ê = Δt*(sum(e .* cis.(v .* s)) - (e[1] + e[end]*cis(v*s[end]))/2)
            M = exp(-v^2/2)/sqrt(2π)
            h = echo_filament(v; α = 1.0, k₁ = 1.0)[1]
            worst = max(worst, abs(M/2 + v*M*ê - h)/abs(h))
        end
        @test worst < 2e-5

        # With the field off the theory is the closed form's small-ε limit, to
        # 1.1e-14, which at ε·k₁τ = 0.1 is within J₁'s curvature of the closed
        # form itself: 1.3e-3.
        t = 10.01:0.02:40.0
        setup = (α = 0.01, ε = 0.01, τ = 10.0, k₁ = 1.0, k₂ = 1.5)
        off = echo_second_order(t; setup..., field = false).echo
        linear = [-im*0.01*(0.5*0.01*(u - 10)/2)*exp(-(0.5*u - 15)^2/2) for u in t]
        closed = [echo_closed_form(u; setup...) for u in t]
        @test maximum(abs, off .- linear) < 1e-12*maximum(abs, linear)
        @test maximum(abs, off .- closed) < 3e-3*maximum(abs, closed)
    end
end
