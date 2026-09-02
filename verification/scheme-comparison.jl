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
    Δx = L/64
    x = collect(Δx:Δx:L)
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
    Δx = L/64
    x = collect(Δx:Δx:L)
    v = collect(-6:0.1:6)
    t = collect(0.0:0.05:40.0)
    f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.5*cos(K*x)))'
    return x, v, f₀, t
end

# `fmax = 1.0` brackets the data, whose peak is 1/√(2π) ≈ 0.399 -- the same
# bounds the harness gives its own PFCNonUniform default.
#
# `Godunov(PiecewiseLinear())` without a limiter is missing on purpose: its
# amplification factor exceeds 1 for every mode, which `test_amplification.jl`
# asserts analytically, so it has nothing to contribute here but a NaN.
schemes() = [
    ("Upwind",                 Upwind()),
    ("LaxWendroff",            LaxWendroff()),
    ("Godunov constant",       Godunov(PiecewiseConstant())),
    ("Godunov VanLeer",        Godunov(PiecewiseLinear(), VanLeer())),
    ("SemiLagrangian linear",  SemiLagrangian(LinearSpline())),
    ("SemiLagrangian cubic",   SemiLagrangian(CubicSpline())),
    ("PFC",                    PFC(fmin = 0.0, fmax = 1.0)),
    ("PFCNonUniform",          nothing),          # the harness default
]

function run_one(scheme, grid = landau_grid)
    x, v, f₀, t = grid()
    sx = scheme
    sv = scheme
    elapsed = @elapsed begin
        r = vlasov_poisson(x, v, f₀, t; scheme_x = sx, scheme_v = sv, invariants = true)
    end
    γ, _ = damping_rate(t, r.ε_e; tmin = 6.0, tmax = 30.0)
    ω, _ = oscillation_frequency(t, r.ε_e; tmin = 6.0, tmax = 30.0)
    mass = maximum(abs, r.mass[1:end-1] .- r.mass[1])/r.mass[1]
    l2 = (r.l2[end-1] - r.l2[1])/r.l2[1]
    fmin = minimum(r.fmin[1:end-1])
    return (; γ, ω, elapsed, mass, l2, fmin)
end

function main()
    println("Landau damping at k = ", K, ", by advection scheme")
    println("grid: 64 x 81, Δt = 0.08, t ≤ 70, Courant 0.81 on the fastest row")
    println("analytic: γ = ", γ_ANALYTIC, ", ω = ", ω_ANALYTIC, "\n")

    # One warm-up so the table is not measuring the compiler.
    run_one(nothing)

    results = Tuple{String,NamedTuple}[]
    for (name, scheme) in schemes()
        result = try
            run_one(scheme)
        catch err
            println(@sprintf("%-24s failed: %s", name, sprint(showerror, err)))
            continue
        end
        push!(results, (name, result))
    end

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
            run_one(scheme, nonlinear_grid)
        catch err
            @printf("%-24s failed: %s\n", name, first(sprint(showerror, err), 60))
            continue
        end
        @printf("%-24s %13.2e %12.1e %11s\n",
                name, r.fmin, r.mass, r.fmin ≥ 0 ? "yes" : "NO")
    end

    println("\nThe two schemes that lead on damping-rate accuracy are the two that")
    println("go negative. That is Godunov's theorem showing up in the physics:")
    println("`LaxWendroff` is second order and not monotone, and the cubic spline")
    println("is third order and not monotone, so both overshoot a steep gradient")
    println("and both undershoot below zero on the other side. At 1% amplitude")
    println("that never surfaces; at 50% it does, and in a run that matters it is")
    println("fatal rather than untidy -- a negative f has no entropy, and the")
    println("diagnostic here threw DomainError on both until the integrand was")
    println("guarded. It is what PFC's limiter exists to prevent, and why the")
    println("harness defaults to it despite not leading the first table.")
end

main()
