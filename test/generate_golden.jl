# Regenerates test/data/golden.txt. Not run by the test suite.
#
#     julia --project=. test/generate_golden.jl
#
# Only run this when a numerical change is *intended*, and say so in the commit.
using Vasilek

include(joinpath(@__DIR__, "golden_cases.jl"))

open(joinpath(@__DIR__, "data", "golden.txt"), "w") do io
    println(io, "# scheme  <64 Float64 values as UInt64 bit patterns>")
    println(io, "# regenerate with: julia --project=. test/generate_golden.jl")
    for (name, scheme) in GOLDEN_CASES
        out = golden_run(scheme)
        print(io, name)
        for v in out
            print(io, " 0x", string(reinterpret(UInt64, v), base = 16, pad = 16))
        end
        println(io)
    end
end
println("wrote test/data/golden.txt")
