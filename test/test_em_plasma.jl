@isdefined(vlasov_poisson) || include(joinpath(@__DIR__, "verification_harness.jl"))

using Vasilek.FDTD1D

"""
The dispersion relation of an electromagnetic wave in uniform plasma, which is
the one thing the transverse side of the reduced model does and the one thing
nothing measured.

`test_fdtd_1d.jl` pins the Yee update in vacuum, to a closed form, mode by mode.
The plasma current that turns that solver into half of a laser-plasma model had
no such test: it was exercised only through `wakefield`, where it enters a
quantity -- the wake -- that is dominated by the *ponderomotive* term. That is
not a small gap. Before the coupling existed at all, the wakefield study's
numbers came out bit-identical whether the transverse current was right, wrong
by a factor of `Δt`, or wrong by thirty-two orders of magnitude; the study can
now see the current, but only through one figure of merit at one set of
parameters.

Here the same `transverse_step!` is driven on its own, in a plasma of constant
density, and the frequency of a standing mode is compared against
[`em_omega`](@ref) -- exact for this update, not the continuum `ωₚ² + k²`. Sign,
`Δt` factor and density all enter that relation, and the last testset runs the
sign the other way to show what a failure looks like.

A standing PEC mode is the right probe for the same reason it is in the vacuum
test: `ey[1]` and `ey[end]` are held at zero by the boundary condition, so
`sin(kx)` with `k = mπ/L` is an exact eigenfunction of the discrete operator,
and with a uniform density the current stays proportional to it. The frequency
is read off the mode's own projection `Σ ey·sin(kx)` rather than a point sample,
which the vacuum test learned the hard way: a probe at `L/4` is exactly on a
node for every `m` divisible by four.
"""

const EMP_NO_PULSE = (y = (t, x) -> 0.0, z = (t, x) -> 0.0)
const EMP_N = 200
const EMP_Δx = 2π/20                      # twenty cells per vacuum wavelength
const EMP_STEPS = 4000

"""
    em_plasma_mode(; density, cfl, m, steps, current_sign)

Seed `sin(kx)` into `ey` at rest and step it with [`transverse_step!`](@ref),
returning the mode projection history and `(k, Δt)`.

`current_sign = -1` is the physical one, the sign `transverse_step!` applies;
`+1` reverses the current and is used only to show the test has teeth.
"""
function em_plasma_mode(; density, cfl, m, steps = EMP_STEPS, current_sign = -1.0)
    Δx, N = EMP_Δx, EMP_N
    Δt = cfl*Δx
    L = N*Δx
    k = m*π/L
    em = FDTD1D.YeeMesh1D{Float64}(N)
    shape = [sin(k*i*Δx) for i = 0:N]
    em.ey .= shape
    advance! = FDTD1D.make_advance_fields(em, cfl, EMP_NO_PULSE, Δt, Δx, 0.0,
                                          FDTD1D.PML(0, 1.0, Δx, Δt))
    n = fill(density, N + 1)
    pʸ, pᶻ = zeros(N + 1), zeros(N + 1)
    projection = Vector{Float64}(undef, steps)
    for s in 1:steps
        if current_sign < 0
            transverse_step!(advance!, em, pʸ, pᶻ, n, s*Δt, Δt)
        else
            # The same three lines with the current the other way round.
            pʸ .= pʸ .+ em.ey.*Δt
            pᶻ .= pᶻ .+ em.ez.*Δt
            advance!(s*Δt, (y = pʸ.*n.*Δt, z = pᶻ.*n.*Δt))
        end
        projection[s] = sum(em.ey .* shape)
    end
    return projection, k, Δt
end

"Frequency from the mean spacing of upward zero crossings, as in the vacuum test."
function crossing_frequency(projection, Δt)
    ups = [s for s in 2:length(projection) if projection[s-1] < 0 ≤ projection[s]]
    length(ups) ≥ 3 || error("crossing_frequency: $(length(ups)) crossings, need 3")
    return 2π/((ups[end] - ups[1])/(length(ups) - 1)*Δt), length(ups)
end

@testset "An EM wave in uniform plasma obeys the discrete dispersion relation" begin
    # Measured relative departure from `em_omega` over every case below:
    # worst 3.2e-4, and it is the estimator rather than the scheme -- a crossing
    # is located only to within Δt, so the shortest runs in wall-clock terms
    # (cfl = 0.05, two dozen crossings) are the least resolved. Held to 1e-3,
    # the same tolerance the vacuum dispersion test uses for the same reason.
    worst = 0.0
    for density in (0.1, 0.3), cfl in (0.05, 0.5, 0.9), m in (1, 5, 20, 50)
        projection, k, Δt = em_plasma_mode(; density = density, cfl = cfl, m = m)
        ω, ncross = crossing_frequency(projection, Δt)
        ω_theory = em_omega(k, density, EMP_Δx, Δt)
        dev = abs(ω - ω_theory)/ω_theory
        worst = max(worst, dev)
        @test dev < 1e-3
        # Above the cutoff, always: an EM wave cannot exist below `√n`, and the
        # discrete relation says so through the `nΔt²/4` term.
        @test ω > sqrt(density)
        if m == 50
            println("  n = ", density, " cfl = ", rpad(cfl, 5), " m = ", lpad(m, 2),
                    ": ω = ", rpad(round(ω; digits = 6), 9),
                    " vs ", rpad(round(ω_theory; digits = 6), 9),
                    " (", round(100*dev; sigdigits = 2), "%, ", ncross, " crossings)")
        end
    end
    println("  worst departure from (2/Δt)²sin²(ωΔt/2) = (2/Δx)²sin²(kΔx/2) + n: ", worst)

    @testset "the plasma cutoff is there" begin
        # At `k → 0` the relation collapses to the discrete plasma frequency,
        # `(2/Δt)·asin(√n·Δt/2)`, which is what a wave at the cutoff does: it
        # oscillates in place. The longest mode this grid holds has
        # `k = 0.05`, so `k²/n` is 0.8% at n = 0.3 and the measured frequency
        # should sit just above `√n`.
        for density in (0.1, 0.3)
            projection, k, Δt = em_plasma_mode(; density = density, cfl = 0.5, m = 1)
            ω, _ = crossing_frequency(projection, Δt)
            println("  n = ", density, ": longest mode (k = ", round(k; digits = 4),
                    ") oscillates at ω = ", round(ω; digits = 5),
                    " against √n = ", round(sqrt(density); digits = 5),
                    " and the discrete ωₚ = ",
                    round(2/Δt*asin(sqrt(density)*Δt/2); digits = 5))
            @test isapprox(ω, sqrt(density + k^2); rtol = 0.01)
        end
    end

    @testset "with no plasma it is the vacuum relation" begin
        # `em_omega` must reduce to `sin(ωΔt/2) = cfl·sin(kΔx/2)`, which
        # `test_fdtd_1d.jl` already holds the solver to. Checking the formula
        # rather than the run, because the run is checked there.
        for cfl in (0.5, 0.9, 1.0), m in (1, 20, 90)
            k = m*π/(EMP_N*EMP_Δx)
            Δt = cfl*EMP_Δx
            @test em_omega(k, 0.0, EMP_Δx, Δt) ≈ 2/Δt*asin(cfl*sin(k*EMP_Δx/2))
        end
    end

    @testset "and the other sign of the current is not a small error" begin
        # `docs/normalization.md` says flipping the sign turns the oscillation
        # into exponential growth. Asserted here rather than described: the same
        # mode with the current reversed grows without bound, so the sign in
        # `transverse_step!` is load-bearing and not a convention one could take
        # either way.
        #
        # Measured after 4000 steps at n = 0.3, cfl = 0.5, m = 5: the projection
        # reaches 4.8e+134 against the 1.0e+2 it starts at and stays at with the
        # physical sign.
        good, _, _ = em_plasma_mode(; density = 0.3, cfl = 0.5, m = 5)
        bad, _, _ = em_plasma_mode(; density = 0.3, cfl = 0.5, m = 5, current_sign = +1.0)
        println("  peak projection: physical sign ", round(maximum(abs, good); sigdigits = 3),
                ", reversed ", round(maximum(abs, bad); sigdigits = 3))
        @test maximum(abs, good) < 2*abs(good[1])
        @test maximum(abs, bad) > 1e10*maximum(abs, good)
    end
end
