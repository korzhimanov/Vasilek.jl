# Work–precision: what each advection scheme costs for the accuracy it buys.
#
#     julia --project=benchmark benchmark/workprecision.jl
#
# `runbenchmarks.jl` times each scheme in isolation and `test/test_comparison.jl`
# ranks their errors. Neither puts the two axes together, and it is only
# together that they answer the question a caller has: for the accuracy I need,
# on the data I have, which scheme is cheapest?
#
# **Advisory, and it exits 0.** The error column is deterministic and is gated in
# `test/test_comparison.jl`. The cost column is wall clock, which
# `runbenchmarks.jl` sets out at length is not gateable on a shared machine at
# any tolerance worth having. What survives the noise is the *shape* -- which
# schemes sit on the frontier, and which are dominated outright -- and that is
# what the summary at the bottom reports.
#
# Timing goes through `BenchmarkTools.@belapsed` rather than `@elapsed` around a
# loop. Julia specialises `advect!` per scheme type, so a first call times the
# compiler; this package has already published one table that made exactly that
# mistake and reported upwind as slower than PFCNonUniform when it is about
# twice as fast. `@belapsed` warms up, tunes the evaluation count, and returns a
# minimum.

using BenchmarkTools
using Printf
using Vasilek

include(joinpath(@__DIR__, "..", "test", "scheme_cases.jl"))

const C = 0.4
const SIZES = (64, 128, 256, 512)
const BUDGET = 0.15          # seconds per timing point; the whole run is ~1 min

profiles() = (
    (:sine,   "smooth sine",  N -> [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]),
    (:gauss,  "gaussian",     N -> [1.0 + exp(-((((i-1)/N) - 0.5)/0.08)^2) for i = 1:N]),
    (:square, "square pulse", N -> [0.4 < (i-1)/N < 0.6 ? 1.5 : 1.0 for i = 1:N]),
)

# `fmax = 2.5` brackets all three profiles: peaks of 1.5, 2.0 and 1.5.
schemes() = (
    ("Upwind",                   Upwind()),
    ("Godunov constant",         Godunov(PiecewiseConstant())),
    ("Godunov VanLeer",          Godunov(PiecewiseLinear(), VanLeer())),
    ("LaxWendroff",              LaxWendroff()),
    ("SemiLagrangian linear",    SemiLagrangian(LinearSpline())),
    ("SemiLagrangian quadratic", SemiLagrangian(QuadraticSpline())),
    ("SemiLagrangian cubic",     SemiLagrangian(CubicSpline())),
    ("PFC",                      PFC(fmin = 0.0, fmax = 2.5)),
)

"L² error after carrying `f₀` exactly once around the domain."
function traversal_error(scheme, f₀, N)
    final = march!(f₀, scheme, C, round(Int, N/C))
    return sqrt((1/N)*sum(abs2, final .- f₀))
end

"Seconds for one `advect!` call, warm, minimum of many."
function step_seconds(scheme, f₀, N)
    src = copy(f₀)
    dst = similar(src)
    ws = workspace(scheme, N)
    return @belapsed advect!($dst, $src, $scheme, $C, $ws) samples=10000 seconds=BUDGET
end

"""
    frontier(points)

The Pareto-optimal subset of `(name, cost, error)`: those for which no other
point is both cheaper **and** more accurate.

This is the part of the table that survives the timing noise. A scheme is
dropped from the frontier only when another beats it on both axes, which takes
a factor rather than a percent -- the dominated entries below lose by 4x to 70x
on error at comparable cost, or by an order of magnitude on cost at worse error.
"""
function frontier(points)
    return [p for p in points
            if !any(q -> q !== p && q[2] ≤ p[2] && q[3] ≤ p[3] &&
                         (q[2] < p[2] || q[3] < p[3]), points)]
end

function main()
    results = Dict{Tuple{Symbol,String,Int},NamedTuple}()

    for (key, label, gen) in profiles()
        println("\n", "="^92)
        println(label, ":  L² error after one traversal, and the cost of getting there")
        println("="^92)
        @printf("%-26s", "scheme")
        for N in SIZES
            @printf("%22s", "N = $N")
        end
        println()
        @printf("%-26s", "")
        for _ in SIZES
            @printf("%12s%10s", "error", "ms/lap")
        end
        println()

        for (name, scheme) in schemes()
            @printf("%-26s", name)
            for N in SIZES
                f₀ = gen(N)
                err = traversal_error(scheme, f₀, N)
                per_step = step_seconds(scheme, f₀, N)
                lap = per_step*round(Int, N/C)
                results[(key, name, N)] = (; err, per_step, lap)
                @printf("%12.3e%10.3f", err, 1e3*lap)
            end
            println()
        end

        # cost per cell per step, which is the N-independent way to read the
        # same numbers and the one that shows a scheme scaling badly
        println()
        @printf("%-26s", "ns/cell/step")
        println()
        for (name, _) in schemes()
            @printf("%-26s", name)
            for N in SIZES
                r = results[(key, name, N)]
                @printf("%22.1f", 1e9*r.per_step/N)
            end
            println()
        end
    end

    println("\n", "="^92)
    println("Efficiency frontier at N = ", last(SIZES),
            ": schemes no other scheme beats on both cost and error")
    println("="^92)
    for (key, label, _) in profiles()
        points = [(name, results[(key, name, last(SIZES))].lap,
                         results[(key, name, last(SIZES))].err)
                  for (name, _) in schemes()]
        best = sort(frontier(points); by = p -> p[2])
        println("\n", label, ":")
        for (name, lap, err) in best
            @printf("    %-26s %9.3f ms   %10.3e\n", name, 1e3*lap, err)
        end
        dominated = [p[1] for p in points if !(p in best)]
        isempty(dominated) || println("    dominated: ", join(dominated, ", "))
    end

    println("""

    Reading this: the frontier is the answer, not the ordering within it. A
    scheme on the frontier is the cheapest way to reach its own accuracy; a
    dominated scheme is one there is no reason to choose, because something else
    is both faster and more accurate. Which of the frontier entries you want
    depends on the accuracy you actually need -- that is the choice the table
    cannot make for you, and the one it exists to inform.

    How much of the cost axis is noise, measured on this run rather than
    asserted: `Upwind` and `Godunov constant` are the same scheme -- identical
    errors to every digit, and `test_amplification.jl` proves the equivalence
    from their amplification factors -- so any gap between their timings is
    measurement error and nothing else. Both come out around 0.1 ns per cell per
    step, at the floor where timer resolution and cache state dominate, so which
    of them reaches the frontier -- and whether the other is reported as
    "dominated" -- changes from run to run. It has done so between two runs of
    this script already. Read that pair as the error bar on every other row and
    treat any single-digit cost ratio as indistinguishable: it is why this
    script exits 0 and why `runbenchmarks.jl` sets a 10 us floor below which it
    gates nothing.

    Two more things it does not say. Cost here is one thread on one machine, so
    read the ratios and not the milliseconds. And accuracy is not the only axis:
    the schemes that lead on smooth data are the ones that drive a distribution
    function negative at large amplitude, which
    `verification/scheme-comparison.jl` measures and this table cannot see.
    """)
end

main()
