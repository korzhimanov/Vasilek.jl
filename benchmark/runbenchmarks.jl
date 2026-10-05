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

# `haskey` on a BenchmarkGroup does not follow a multi-key path, so presence is
# tested by indexing.
has_leaf(group, key) = try group[key]; true catch; false end

if isfile(PARAMS_FILE)
    stored = BenchmarkTools.load(PARAMS_FILE)[1]
    loadparams!(SUITE, stored, :evals, :samples)
    # `loadparams!` skips entries the file does not have, which then ran at the
    # default evals = 1 -- how a new entry went untuned. Tune those here and
    # say so.
    untuned = [k for (k, _) in BenchmarkTools.leaves(SUITE) if !has_leaf(stored, k)]
    for k in untuned
        tune!(SUITE[k])
    end
    isempty(untuned) || println("tuned ", length(untuned), " entries missing from ",
                                PARAMS_FILE, ": ", join(join.(untuned, " / "), ", "),
                                "\nre-run with --rebaseline to store them")
    REBASELINE && BenchmarkTools.save(PARAMS_FILE, params(SUITE))
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

"""
    baseline_versions(file)

The `Julia` and `BenchmarkTools` versions `BenchmarkTools.save` wrote at the head
of `file`, which `BenchmarkTools.load` reads past. Matched rather than parsed:
the environment has no JSON package of its own.
"""
function baseline_versions(file)
    text = read(file, String)
    version(key) = (m = match(Regex("\"$key\":\"([^\"]*)\""), text)) === nothing ? "unknown" : m[1]
    return version("Julia"), version("BenchmarkTools")
end

# The stored baseline was measured on the development machine, and
# `results.json` records the Julia and BenchmarkTools it ran under in its header,
# which `--rebaseline` rewrites. On any other machine or Julia the comparison
# below measures that difference as much as any change in the code, so the
# header is printed beside this run's.
baseline = BenchmarkTools.load(RESULTS_FILE)[1]
let (julia, tools) = baseline_versions(RESULTS_FILE)
    println("\nbaseline ", RESULTS_FILE, ": Julia ", julia, ", BenchmarkTools ", tools,
            "; this run: Julia ", VERSION, ", BenchmarkTools ", pkgversion(BenchmarkTools))
end
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
    unbased = String[]
    for (key, trial) in BenchmarkTools.leaves(results)
        if !haskey(basetimes, key)
            push!(unbased, join(key, " / "))
            continue
        end
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
    return offenders, ignored, unbased
end

offenders, ignored, unbased = gate(results, baseline)

println("\n", ignored, " of ", length(BenchmarkTools.leaves(results)),
        " entries are below the ", round(Int, GATE_FLOOR_NS/1000),
        " us floor: reported above, not gated.")
isempty(unbased) || println(length(unbased), " entries have no baseline and were not compared: ",
                            join(unbased, ", "), "\nre-run with --rebaseline to record them.")

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
