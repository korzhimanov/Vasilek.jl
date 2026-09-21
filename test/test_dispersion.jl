@isdefined(Z) || include(joinpath(@__DIR__, "dispersion.jl"))

"""
The kinetic dispersion relation, checked before anything is allowed to be
measured against it.

`dispersion.jl` is a reference implementation: the physics tests compare fitted
rates against the roots it returns, so an error here would move a target rather
than fail a test, and every number in the suite that leans on it would drift in
the same direction at once. That is the one failure mode a tolerance cannot
catch, which is why the checks below are identities and independent quadratures
rather than agreement with stored values.

**The exception is `LANDAU_ROOTS`**, which is exactly a stored value -- the
three roots this suite carried as constants before there was a solver. Keeping
them here, as a fixture for the solver rather than as the target of a physics
test, is what makes the replacement auditable: the computed roots have to
reproduce the numbers every Landau assertion was held to until now.

Cheap enough for the default test run: the whole file is a few hundred `erfcx`
calls.
"""

"""
Roots of the Landau dispersion relation for a Maxwellian, `(γ, ω_r)` by `kλ_D`.

The tabulated values `test_verification.jl` carried before `landau_root` existed.
They are kept as the fixture the solver is checked against, and are not used as
a target anywhere else.
"""
const LANDAU_ROOTS = Dict(0.3 => (0.01262, 1.15985),
                          0.4 => (0.06613, 1.28506),
                          0.5 => (0.15336, 1.41566))

"Trapezoid over a window wide enough that the Gaussian tail is below round-off."
function quad(f; lo = -12.0, hi = 12.0, n = 24001)
    h = (hi - lo)/(n - 1)
    s = 0.5*(f(lo) + f(hi))
    for i in 2:n-1
        s += f(lo + (i-1)*h)
    end
    return s*h
end

@testset "The plasma dispersion function" begin
    @testset "Z at zero and on the real axis" begin
        # Z(0) = i√π exactly, and the imaginary part on the real axis is
        # √π·exp(-ζ²) -- the Landau term, 2.8e-28 at ζ = 8. The unscaled
        # `exp(-ζ²)erfc(-iζ)` still gets that right: the two forms agree to the
        # last bit up to ζ = 26. `erfcx` earns its place past that, where
        # `erfc(-iζ)` overflows and the unscaled form returns `Re Z = -Inf` --
        # which is why ζ = 30 is in the asymptotic loop below.
        @test Z(0.0 + 0.0im) ≈ im*sqrt(π) rtol = 1e-15
        for ζ in (0.5, 2.0, 5.0, 8.0)
            @test imag(Z(complex(ζ, 0.0))) ≈ sqrt(π)*exp(-ζ^2) rtol = 1e-12
        end

        # and the asymptotic expansion of the real part, which pins the branch:
        # Re Z → -1/ζ - 1/(2ζ³) - 3/(4ζ⁵), the next term being -15/(8ζ⁷).
        for ζ in (5.0, 8.0, 30.0)
            series = -1/ζ - 1/(2ζ^3) - 3/(4ζ^5)
            @test isfinite(Z(complex(ζ, 0.0)))
            @test isapprox(real(Z(complex(ζ, 0.0))), series; atol = 2*15/(8ζ^7))
        end
    end

    @testset "Z is the Cauchy integral where the integral converges" begin
        # The definition -- (1/√π)∫exp(-t²)/(t - ζ)dt -- holds only for
        # Im ζ > 0, which is why `Z` is written through the entire function `w`
        # instead. Above the axis the two must agree, and that is what says the
        # continuation below the axis, where every damped root lives, starts
        # from the right function.
        #
        # Measured worst |Δ| over the four points below: 3.5e-15. The quadrature
        # is a plain trapezoid, which is spectrally accurate here -- the
        # integrand is analytic and Gaussian-damped, and the nearest pole sits
        # half a unit off the contour.
        worst = 0.0
        for ζ in (0.5 + 0.5im, -0.3 + 0.8im, 2.0 + 1.0im, 0.0 + 2.0im)
            cauchy = quad(t -> exp(-t^2)/(t - ζ))/sqrt(π)
            worst = max(worst, abs(Z(ζ) - cauchy))
            @test Z(ζ) ≈ cauchy rtol = 1e-10
        end
        println("  Z vs its Cauchy integral: worst |Δ| = ", worst)
    end

    @testset "Zprime is the identity it claims" begin
        # Z'(ζ) = -2(1 + ζZ(ζ)) against a central difference, which is the one
        # line in `dispersion.jl` a typo would leave looking reasonable.
        #
        # Measured worst relative departure: 9.3e-11, which is the difference's
        # own truncation error (h²|Z'''|/6 at h = 1e-5), not the identity's.
        # Hence `rtol = 1e-8` rather than machine precision.
        worst = 0.0
        h = 1e-5
        for ζ in (0.3 + 0.2im, 1.5 - 0.4im, 2.8 - 0.1im, -1.0 + 0.6im)
            fd = (Z(ζ + h) - Z(ζ - h))/(2h)
            worst = max(worst, abs(Zprime(ζ) - fd)/abs(fd))
            @test isapprox(Zprime(ζ), fd; rtol = 1e-8)
        end
        println("  Z' vs central difference: worst relative |Δ| = ", worst)
    end

    @testset "the susceptibility is the integral it stands for" begin
        # χ = -(1/k²)∫M'(v)/(v - ω/k)dv for a Maxwellian M, which integrates by
        # parts into (1/k²)(1 + ζZ(ζ)). The quadrature shares no code with the
        # closed form: it never calls `Z`. Measured worst relative departure
        # over the three (ω, k) below: 1.3e-14.
        M(v) = exp(-0.5*v^2)/sqrt(2π)
        M′(v) = -v*M(v)
        worst = 0.0
        for (ω, k) in ((1.4 + 0.3im, 0.5), (1.2 + 0.8im, 0.4), (0.9 + 1.5im, 0.8))
            numeric = -quad(v -> M′(v)/(v - ω/k))/k^2
            worst = max(worst, abs(susceptibility(ω, k) - numeric)/abs(numeric))
            @test isapprox(susceptibility(ω, k), numeric; rtol = 1e-10)
        end
        println("  χ vs quadrature: worst relative |Δ| = ", worst)

        # Two identical half-density species are one full one: this is the only
        # statement about `dielectric`'s sum, and the two-stream root below is
        # built entirely on it.
        one = (density = 1.0, drift = 0.0, vt = 1.0)
        halves = ((density = 0.5, drift = 0.0, vt = 1.0), (density = 0.5, drift = 0.0, vt = 1.0))
        @test dielectric(1.4 + 0.3im, 0.5, (one,)) ≈ dielectric(1.4 + 0.3im, 0.5, halves)
    end
end

@testset "Landau roots reproduce the table they replace" begin
    # The three constants every Landau assertion in this suite was held to,
    # recomputed. Agreement to every digit the table quotes is the point: the
    # replacement has to be a strictly better source for the same numbers, not
    # a different set of numbers that also happens to pass.
    for k in sort(collect(keys(LANDAU_ROOTS)))
        γa, ωa = LANDAU_ROOTS[k]
        ω = landau_root(k)
        println("  k = ", k, "   computed ", round(real(ω); digits = 6), " - ",
                round(-imag(ω); digits = 6), "im   table ", ωa, ", ", γa)
        @test isapprox(real(ω), ωa; rtol = 1e-5)
        @test isapprox(-imag(ω), γa; rtol = 1e-4)
        @test imag(ω) < 0                      # damped, not growing
    end

    # And the root solves the relation it came from, which is what says the
    # secant converged to a root rather than to wherever it ran out of steps.
    for k in (0.3, 0.4, 0.5, 0.7, 1.0)
        residual = abs(dielectric(landau_root(k), k, ((density = 1.0, drift = 0.0, vt = 1.0),)))
        @test residual < 1e-10
    end
end

@testset "Warm two-stream roots" begin
    @testset "the cold limit is the cold closed form" begin
        # `γ_cold` is the vt → 0 limit of this function, so the two must meet
        # there. They do not meet at the vt = 0.3 the runs use, which is the
        # whole reason the warm root exists. Measured at vt = 0.02: -0.006%,
        # -0.005% and +0.002% at a = 0.4, 0.6 and 0.8 -- the residue is the
        # remaining temperature, not the solver.
        γ_cold(a) = sqrt(max(0.0, -((2a^2 + 1) - sqrt(8a^2 + 1))/2))
        for a in (0.4, 0.6, 0.8)
            warm = two_stream_warm(a; vt = 0.02)
            println("  a = ", a, "  vt = 0.02: warm ", round(warm; digits = 5),
                    " vs cold ", round(γ_cold(a); digits = 5),
                    "  (", round(100*(warm - γ_cold(a))/γ_cold(a); digits = 3), "%)")
            @test isapprox(warm, γ_cold(a); rtol = 0.01)
        end

        # The departure grows with temperature, monotonically, at fixed a.
        rates = [two_stream_warm(0.6; vt = vt) for vt in (0.02, 0.2, 0.3, 0.45, 0.6)]
        @test issorted(rates; rev = true)
    end

    @testset "warm beams outgrow cold ones near the band edge" begin
        # The cold branch is identically zero for a ≥ 1 while warm beams are
        # still unstable there, so by continuity the warm rate must cross above
        # the cold one somewhere below the edge. Measured, the crossing sits
        # between a = 0.75 and a = 0.8 at vt = 0.3.
        γ_cold(a) = sqrt(max(0.0, -((2a^2 + 1) - sqrt(8a^2 + 1))/2))
        @test two_stream_warm(0.6) < γ_cold(0.6)
        @test two_stream_warm(0.75) < γ_cold(0.75)
        @test two_stream_warm(0.8) > γ_cold(0.8)
        @test two_stream_warm(0.9) > γ_cold(0.9)
        for a in (0.6, 0.75, 0.8, 0.9, 1.0)
            println("  a = ", rpad(a, 4), " warm ", rpad(round(two_stream_warm(a); digits = 5), 7),
                    " cold ", round(γ_cold(a); digits = 5))
        end

        # At the cold boundary the warm branch is what the suite's `a = 1.0`
        # case grows at, and beyond it the cold form has nothing to say at all.
        @test two_stream_warm(1.0) > 0.05
    end

    @testset "and stop somewhere" begin
        # Bisection reports 0.0 rather than a number when there is no growing
        # root in (0, 1]. Cold-stable a = 1.2 is warm-unstable or not depending
        # on vt, and the assertion here is the shape of that: cold enough beams
        # at a fixed a beyond the boundary are stable.
        @test two_stream_warm(1.2; vt = 0.05) == 0.0
        @test two_stream_warm(2.0; vt = 0.3) == 0.0

        # A returned rate is always a root, never a bisection artefact.
        for (a, vt) in ((0.4, 0.3), (0.8, 0.3), (1.0, 0.3), (0.6, 0.6))
            γ = two_stream_warm(a; vt = vt)
            beams = ((density = 0.5, drift = 3.0, vt = vt),
                     (density = 0.5, drift = -3.0, vt = vt))
            @test abs(dielectric(im*γ, a/3.0, beams)) < 1e-6
        end
    end
end
