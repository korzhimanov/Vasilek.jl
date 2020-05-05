using BenchmarkTools

const SUITE = BenchmarkGroup()
const PARAMS_FILE = joinpath(dirname(@__FILE__), "params.json")
const RESULTS_FILE = joinpath(dirname(@__FILE__), "results.json")

include(joinpath(dirname(@__FILE__), "MaxwellBenchmarks.jl"))
using .MaxwellBenchmarks
SUITE["maxwell"] = MaxwellBenchmarks.SUITE

if isfile(PARAMS_FILE)
    loadparams!(SUITE, BenchmarkTools.load(PARAMS_FILE)[1])
else
    tune!(SUITE)
    BenchmarkTools.save(PARAMS_FILE, params(SUITE))
end

if isfile(RESULTS_FILE)
    m = BenchmarkTools.load(RESULTS_FILE)[1]
    r = run(SUITE)
    j = judge(minimum(m), minimum(r))
    println(j)
    if !isregression(j)
        println("Saving new results")
        BenchmarkTools.save(RESULTS_FILE, r)
    end
else
    r = run(SUITE)
    print(minimum(r))
    BenchmarkTools.save(RESULTS_FILE, r)
end
