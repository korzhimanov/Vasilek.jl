using BenchmarkTools

const SUITE = BenchmarkGroup()
const PARAMS_FILE = joinpath(dirname(@__FILE__), "params.json")
const RESULTS_FILE = joinpath(dirname(@__FILE__), "results.json")

include(joinpath(dirname(@__FILE__), "MaxwellBenchmarks.jl"))
using .MaxwellBenchmarks
SUITE["Maxwell"] = MaxwellBenchmarks.SUITE

include(joinpath(dirname(@__FILE__), "VlasovBenchmarks.jl"))
using .VlasovBenchmarks
SUITE["Vlasov"] = VlasovBenchmarks.SUITE

if isfile(PARAMS_FILE)
    loadparams!(SUITE, BenchmarkTools.load(PARAMS_FILE)[1])
else
    tune!(SUITE)
    BenchmarkTools.save(PARAMS_FILE, params(SUITE))
end

r = run(SUITE["Vlasov"])
println(minimum(r))

if isfile(RESULTS_FILE)
    m = BenchmarkTools.load(RESULTS_FILE)[1]
    j = judge(minimum(r), minimum(m))
    println(j)
    if !isregression(j) && isimprovement(j)
        println("Saving new results")
        BenchmarkTools.save(RESULTS_FILE, r)
    end
else
    BenchmarkTools.save(RESULTS_FILE, r)
end
