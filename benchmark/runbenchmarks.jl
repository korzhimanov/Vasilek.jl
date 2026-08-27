# Benchmark driver.
#
#     julia --project=benchmark benchmark/runbenchmarks.jl
#     julia --project=benchmark benchmark/runbenchmarks.jl --strict
#     julia --project=benchmark benchmark/runbenchmarks.jl --rebaseline
#
# Advisory by default: it reports how the suite compares with the stored
# baseline and exits 0. --strict turns a regression into a non-zero exit.
# See the note on tolerances below for why that is not the default, and why
# this is not wired into per-PR CI.
using BenchmarkTools

const SUITE = BenchmarkGroup()
const PARAMS_FILE = joinpath(@__DIR__, "params.json")
const RESULTS_FILE = joinpath(@__DIR__, "results.json")
const REBASELINE = "--rebaseline" in ARGS
const STRICT = "--strict" in ARGS

# Measured, not guessed. Repeated runs of the whole suite against a freshly
# written baseline, on an idle development machine with nothing changed:
#
#   median deviation 2.2%, p90 about 6%
#   20-50% on sub-microsecond entries (timer resolution, cache state)
#   26% on a 13 us entry, 25-32% on 1.5 us entries
#
# At BenchmarkTools default 5% an unchanged re-run reports regressions. So
# does 25%. Adding a 1 us floor did not help; a 10 us floor still tripped
# once in four runs. Wall-clock gating of kernels this fast is not reliable
# on a machine that is not isolated and pinned, at any tolerance worth
# having.
#
# So: the comparison is advisory by default and exits 0. --strict makes it
# fail, for someone running on a quiet machine who wants that. The
# deterministic gate that does block CI is the allocation test.
const TIME_TOLERANCE = 0.15
const GATE_FLOOR_NS = 10_000.0

include(joinpath(@__DIR__, "MaxwellBenchmarks.jl"))
using .MaxwellBenchmarks
SUITE["Maxwell"] = MaxwellBenchmarks.SUITE

include(joinpath(@__DIR__, "VlasovBenchmarks.jl"))
using .VlasovBenchmarks
SUITE["Vlasov"] = VlasovBenchmarks.SUITE

if isfile(PARAMS_FILE)
    loadparams!(SUITE, BenchmarkTools.load(PARAMS_FILE)[1], :evals, :samples)
else
    tune!(SUITE)
    BenchmarkTools.save(PARAMS_FILE, params(SUITE))
end

# The whole suite. This used to run SUITE["Vlasov"] only, so the Maxwell group
# was tuned and then never executed.
results = minimum(run(SUITE))
println(results)

if REBASELINE || !isfile(RESULTS_FILE)
    # Saving only on an improvement lets the baseline ratchet downwards and
    # never back, so an intentional slowdown could not be recorded. Saving is
    # now an explicit request. `minimum(...)` keeps the estimate rather than
    # every sample: the file was 4.4 MB of raw trials, and is now about 13 kB.
    BenchmarkTools.save(RESULTS_FILE, results)
    println("\nbaseline written to ", RESULTS_FILE)
    exit(0)
end

baseline = BenchmarkTools.load(RESULTS_FILE)[1]
println("\n", judge(results, baseline; time_tolerance = TIME_TOLERANCE))

"""
    gate(results, baseline)

Entries slower than the floor and worse than the tolerance, plus a count of
those ignored for being too fast to measure. A function rather than top-level
code: accumulating into a bare loop variable at script scope hits Julia's soft
scope rule and fails with `UndefVarError`.
"""
function gate(results, baseline)
    basetimes = Dict(k => BenchmarkTools.time(v) for (k, v) in BenchmarkTools.leaves(baseline))
    offenders = Tuple{String,Float64,Float64}[]
    ignored = 0
    for (key, trial) in BenchmarkTools.leaves(results)
        haskey(basetimes, key) || continue
        before = basetimes[key]
        if before < GATE_FLOOR_NS
            ignored += 1
            continue
        end
        now = BenchmarkTools.time(trial)
        if now > before*(1 + TIME_TOLERANCE)
            push!(offenders, (join(key, " / "), before, now))
        end
    end
    return offenders, ignored
end

offenders, ignored = gate(results, baseline)

println("\n", ignored, " of ", ignored + length(BenchmarkTools.leaves(results)) - ignored,
        " entries are below the ", round(Int, GATE_FLOOR_NS/1000),
        " us floor: reported above, not gated.")

if isempty(offenders)
    println("no regression above the floor")
    exit(0)
end

println("\nSlower than the stored baseline by more than ",
        round(Int, 100*TIME_TOLERANCE), "%:")
for (name, before, now) in offenders
    println("  ", rpad(name, 46), round(before/1000; digits = 2), " us -> ",
            round(now/1000; digits = 2), " us  (+",
            round(100*(now/before - 1); digits = 1), "%)")
end
println("\nIf it is intended, re-run with --rebaseline.")

if STRICT
    exit(1)
end
println("Advisory only. Pass --strict to make this fail the run.")
