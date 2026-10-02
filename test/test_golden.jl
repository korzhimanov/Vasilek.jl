using Vasilek

include(joinpath(@__DIR__, "golden_cases.jl"))

"""
Bit-for-bit regression net for the advection kernels.

Convergence rates are far too loose to catch an index swapped by one -- that
is a change of a few percent, well inside any order estimate. Bit-identity
is not. This is the safety net the closure-factory refactor and the eventual
scheme-struct rewrite are meant to run against.

Deliberately advection-only. The kernels here are pure `+ - * /`, which IEEE
pins down, and were verified identical across Julia 1.10.12 and 1.12.7 --
except the quadratic and cubic `SemiLagrangian`, whose prefilter is
Interpolations' periodic Woodbury solve. Its operation order belongs to that
package, which the repository does not pin (no Manifest, compat "0.16"), so an
upstream patch could move those two by an ulp with nothing changed here. They
are held to 1e-13 relative instead: still four orders below a one-cell index
error on these data, and blind only to rounding. The
collision operators go through `exp`/`sqrt` and hence libm: BGK already
differs in the last digits between those two versions
(7.83796228241873e-5 vs 7.837962283090717e-5), so a stored bit pattern for
it would be a platform assertion rather than a correctness one.
"""
# Cases whose last bits belong to a dependency's operation order, not to this package.
const GOLDEN_UPSTREAM_ROUNDING = ("SemiLagrangian_quadratic", "SemiLagrangian_cubic")

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
        if name in GOLDEN_UPSTREAM_ROUNDING
            @test isapprox(out, expected; rtol = 1e-13)
        else
            @test out == expected
        end
    end
end
