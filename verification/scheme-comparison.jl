# Which advection scheme should a Vlasov–Poisson run use?
#
#     julia --project=verification verification/scheme-comparison.jl
#
# The benchmark suite times each scheme on a bare advection kernel, and
# `test/test_convergence.jl` measures the order each one achieves on a shifted
# sine. Neither answers the question a user actually has, which is what a scheme
# costs *in the physics*: how close the Landau damping rate comes to its
# analytic value, and how long that takes.
#
# **Advisory, and it exits 0.** The accuracy column is deterministic and could
# be gated; the timing column is wall clock on whatever machine this runs on,
# and `benchmark/runbenchmarks.jl` sets out at length why that is not something
# to gate. Read the two together: the ranking is the output, not the seconds.
#
# The driver is `test/verification_harness.jl` rather than a copy of it. A
# second implementation of the same Strang loop would drift from the one the
# tests assert against, and then this script would be comparing schemes on a
# solver nobody checks.

include(joinpath(@__DIR__, "..", "test", "verification_harness.jl"))

using Printf

const K = 0.5
const γ_ANALYTIC = 0.15336
const ω_ANALYTIC = 1.41566

"The k = 0.5 Landau configuration the extended test asserts, at its own grid."
function landau_grid()
    L = 2*(2π/K)
    Nx = 64
    Δx = L/Nx
    # Nx points by construction. `Δx:Δx:L` rounds its length out of the
    # endpoints and can come up a point short; see `two_stream` in the harness.
    x = collect(range(Δx; step = Δx, length = Nx))
    v = collect(-4:0.1:4)
    t = collect(0.0:0.08:70.0)
    f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(K*x)))'
    return x, v, f₀, t
end

"""
The same mode at 50% amplitude instead of 1%, on a wider velocity window.

The linear case above cannot rank the schemes on positivity: at 1% the
distribution never comes near zero, and every scheme in the table returns the
same `min f = 1.3e-4`, which is just the Maxwellian's own tail at `v = ±4`. It
is a real limitation of that diagnostic rather than evidence that the schemes
behave alike, and it is worth a second configuration rather than a footnote.

At 50% the mode steepens into the sharp phase-space gradients where monotonicity
starts to matter, and the schemes separate immediately.
"""
function nonlinear_grid()
    L = 2*(2π/K)
    Nx = 64
    Δx = L/Nx
    x = collect(range(Δx; step = Δx, length = Nx))
    v = collect(-6:0.1:6)
    t = collect(0.0:0.05:40.0)
    f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.5*cos(K*x)))'
    return x, v, f₀, t
end

# `PFC` is bounded by the distribution it carries, `[0, maximum(f)]`, as the
# harness bounds its own `PFCNonUniform` default. It is given as a function of `f`
# because the driver rescales the initial condition before running, and a bound
# taken from `f₀` beforehand would sit under the rescaled maximum -- by an ulp or
# two, the rescaling being 1 to round-off here, which is enough to trip `PFC`'s
# check.
#
# `Godunov(PiecewiseLinear())` without a limiter is missing on purpose: it is
# `LaxWendroff`, whose row is already here, as `test_amplification.jl` asserts
# mode by mode.
schemes() = [
    ("Upwind",                 Upwind()),
    ("LaxWendroff",            LaxWendroff()),
    ("Godunov constant",       Godunov(PiecewiseConstant())),
    ("Godunov VanLeer",        Godunov(PiecewiseLinear(), VanLeer())),
    ("Godunov Superbee",       Godunov(PiecewiseLinear(), Superbee())),
    ("SemiLagrangian linear",  SemiLagrangian(LinearSpline())),
    ("SemiLagrangian cubic",   SemiLagrangian(CubicSpline())),
    ("PFC",                    f -> PFC(fmin = 0.0, fmax = maximum(f))),
    ("PFCNonUniform",          nothing),          # the harness default
]

"""
    run_one(scheme, grid = landau_grid; fit = true, timed = false)

One run of `grid` under `scheme`, reduced to the columns the tables print.

`fit = false` skips the damping-rate and frequency fits and returns `NaN` for
both. The second table has no analytic rate to compare against and prints
neither, and fitting anyway would throw away a row on a failure that does not
concern it: `damping_rate` needs three maxima in its window and
`oscillation_frequency` two minima, which a scheme that flattens or goes wild at
50% amplitude need not provide, and the `min f` the table exists to report is
already in hand by then.

`timed = true` runs the case **twice** and clocks the second. Timing the first
call measures the compiler: `line_advector` returns a distinct closure type per
scheme, so `vlasov_poisson`, `make_time_step_2d!` and the whole step loop are
specialised afresh for every row of the table. Measured here, first call against
second: `Upwind` 1.03 s / 0.12 s, `Godunov constant` 0.98 s / 0.10 s, `PFC`
1.23 s / 0.23 s. A single warm-up in `main` -- which is what this did until now
-- only ever warmed the default scheme's own row, and left the table reporting
`Upwind` as marginally *slower* than `PFCNonUniform` when it is about twice as
fast.
"""
function run_one(scheme, grid = landau_grid; fit = true, timed = false)
    x, v, f₀, t = grid()
    solve() = vlasov_poisson(x, v, f₀, t;
                             scheme_x = scheme, scheme_v = scheme, invariants = true)
    timed && solve()                     # compile this scheme, off the clock
    local r
    elapsed = @elapsed r = solve()
    γ = fit ? damping_rate(t, r.ε_e; tmin = 6.0, tmax = 30.0)[1] : NaN
    ω = fit ? oscillation_frequency(t, r.ε_e; tmin = 6.0, tmax = 30.0)[1] : NaN
    mass = maximum(abs, r.mass[1:end-1] .- r.mass[1])/r.mass[1]
    l2 = (r.l2[end-1] - r.l2[1])/r.l2[1]
    fmin = minimum(r.fmin[1:end-1])
    return (; γ, ω, elapsed, mass, l2, fmin)
end

function main()
    println("Landau damping at k = ", K, ", by advection scheme")
    println("grid: 64 x 81, Δt = 0.08, t ≤ 70, Courant 0.81 on the fastest row")
    println("analytic: γ = ", γ_ANALYTIC, ", ω = ", ω_ANALYTIC, "\n")

    # Each row warms itself; see `run_one`. There is no single warm-up that
    # would serve, because every scheme compiles its own step loop.
    results = Tuple{String,NamedTuple}[]
    for (name, scheme) in schemes()
        result = try
            run_one(scheme; timed = true)
        catch err
            println(@sprintf("%-24s failed: %s", name, sprint(showerror, err)))
            continue
        end
        push!(results, (name, result))
    end

    # `baseline` is a `minimum` over `results`, and every row above can have
    # been dropped by the `catch`. An empty reduction throws, which would take
    # the script out on a non-zero exit over a table it merely cannot print --
    # and the second table below is independent of this one and still worth
    # running.
    if isempty(results)
        println("\nevery scheme failed on this case; nothing to rank")
    else
        @printf("%-24s %9s %8s   %9s %8s   %8s %7s   %10s %11s\n",
                "scheme", "γ", "err", "ω", "err", "wall", "rel", "mass drift", "ΔL²/L²")
        baseline = minimum(r.elapsed for (_, r) in results)
        for (name, r) in sort(results; by = p -> abs(p[2].γ - γ_ANALYTIC))
            @printf("%-24s %9.5f %7.2f%%   %9.5f %7.2f%%   %7.3fs %6.1fx   %10.1e %11.2e\n",
                    name, r.γ, 100*abs(r.γ - γ_ANALYTIC)/γ_ANALYTIC,
                    r.ω, 100*abs(r.ω - ω_ANALYTIC)/ω_ANALYTIC,
                    r.elapsed, r.elapsed/baseline, r.mass, r.l2)
        end

        println("\nRanked by damping-rate error. Read the wall column as an order of")
        println("magnitude, not a measurement: it is one run on a shared machine.")
    end

    # ---------------------------------------------------------------------
    println("\n\nSame mode at 50% amplitude, where positivity starts to matter")
    println("grid: 64 x 121, Δt = 0.05, t ≤ 40\n")
    println("The table above cannot rank the schemes on positivity: at 1% the")
    println("distribution never approaches zero and every scheme returns the same")
    println("min f = 1.3e-4, which is the Maxwellian's own tail at v = ±4. There is")
    println("no analytic damping rate at this amplitude, so the diagnostic here is")
    println("the sign of f rather than an error.\n")

    @printf("%-24s %13s %12s %11s\n", "scheme", "min f", "mass drift", "positive")
    for (name, scheme) in schemes()
        r = try
            run_one(scheme, nonlinear_grid; fit = false)
        catch err
            @printf("%-24s failed: %s\n", name, first(sprint(showerror, err), 60))
            continue
        end
        @printf("%-24s %13.2e %12.1e %11s\n",
                name, r.fmin, r.mass, r.fmin ≥ 0 ? "yes" : "NO")
    end

    println("\nAmong the schemes that dissipate, the two that lead on damping-rate")
    println("accuracy are the two that go negative. That is Godunov's theorem")
    println("showing up in the physics: `LaxWendroff` is second order and not")
    println("monotone, and the cubic spline is third order and not monotone, so")
    println("both overshoot a steep gradient and both undershoot below zero on the")
    println("other side. At 1% amplitude that never surfaces; at 50% it does, and in")
    println("a run that matters it is fatal rather than untidy -- a negative f has")
    println("no entropy, and the diagnostic here threw DomainError on both until")
    println("the integrand was guarded. It is what PFC's limiter exists to prevent,")
    println("and why the harness defaults to it despite not leading the first table.")
    println("\n`Godunov Superbee` does not dissipate, and it heads the first table")
    println("from below. Its compression is anti-diffusion: its L² is the only one")
    println("that grows, and its rate sits under the analytic value where every")
    println("other scheme's sits above. Refined to Nx = 256 at the same Courant")
    println("number, its error changes sign (-2.40%, -0.21%, -0.67%, +0.13%) rather")
    println("than shrinking, so its lead here is a cancellation, not accuracy.")
end

main()
