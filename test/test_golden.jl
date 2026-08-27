using Vasilek

include(joinpath(@__DIR__, "golden_cases.jl"))

"""
Bit-for-bit regression net for the advection kernels.

Convergence rates are far too loose to catch an index swapped by one -- that
is a change of a few percent, well inside any order estimate. Bit-identity
is not. This is the safety net the closure-factory refactor and the eventual
scheme-struct rewrite are meant to run against.

Deliberately advection-only. The kernels here are pure `+ - * /`, which IEEE
pins down, and were verified identical across Julia 1.10.12 and 1.12.7. The
collision operators go through `exp`/`sqrt` and hence libm: BGK already
differs in the last digits between those two versions
(7.83796228241873e-5 vs 7.837962283090717e-5), so a stored bit pattern for
it would be a platform assertion rather than a correctness one.
"""
function read_golden()
    golden = Dict{String,Vector{Float64}}()
    for line in eachline(joinpath(@__DIR__, "data", "golden.txt"))
        (isempty(line) || startswith(line, "#")) && continue
        parts = split(line)
        golden[parts[1]] = [reinterpret(Float64, parse(UInt64, p)) for p in parts[2:end]]
    end
    return golden
end

@testset "Golden values" begin
    golden = read_golden()
    @test length(golden) == length(GOLDEN_CASES)

    for (name, scheme) in GOLDEN_CASES
        @test haskey(golden, name)
        haskey(golden, name) || continue
        out = golden_run(scheme)
        expected = golden[name]
        if out != expected
            println("GOLDEN MISMATCH $name: max|Δ| = ", maximum(abs, out .- expected),
                    ", max ULP-ish rel = ",
                    maximum(abs.(out .- expected) ./ max.(abs.(expected), eps())))
        end
        @test out == expected
    end
end
