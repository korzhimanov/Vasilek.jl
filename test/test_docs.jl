using Vasilek

# The code in docs/*.md, executed, as `test_readme.jl` does for the README. The
# migration guide's collision example called `workspace(BGK(τ), n)` through a
# name that did not reach the collision operators, and nothing ran it.
#
# A block that shows the 0.1 call and then the 0.2 one is cut at "# 0.2": the
# 0.1 half documents what no longer exists. The blocks use the names a reader's
# code would have, which the preamble defines.
const DOCS_PREAMBLE = """
v  = collect(range(-4, 4; length = 64))
f₀ = exp.(-v.^2 ./ 2) ./ sqrt(2π)
f  = similar(f₀)
c  = 0.4
nsteps = 3
Δx = 0.1
Δt = 0.05
τ  = 1.0
ρ  = sin.(2π .* (0:63) ./ 64)
e  = similar(ρ)
"""

@testset "The docs' examples run" begin
    for file in filter(endswith(".md"), readdir(joinpath(@__DIR__, "..", "docs"); join = true))
        text = read(file, String)
        blocks = [m.captures[1] for m in eachmatch(r"```julia\r?\n(.*?)```"s, text)]
        for code in blocks
            i = findfirst("# 0.2", code)
            current = i === nothing ? code : code[first(i):end]
            occursin("# 0.1", current) && continue      # 0.1 only, nothing to run
            @test (include_string(Module(:DocsExample), "using Vasilek\n" * DOCS_PREAMBLE * current); true)
        end
        println("  ", basename(file), ": ", length(blocks), " julia blocks")
    end
end
