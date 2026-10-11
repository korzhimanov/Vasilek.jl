using Vasilek

# The README's usage examples, executed, so that they cannot drift from the
# API. There are two; a broken fence that hid one would leave fewer.
@testset "The README examples run" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    blocks = [m.captures[1] for m in eachmatch(r"```julia\r?\n(.*?)```"s, readme)]
    println("  julia blocks found in README.md: ", length(blocks))
    @test length(blocks) >= 2
    for code in blocks
        # A fresh module per block, so one cannot lean on names another defined.
        @test (include_string(Module(:READMEExample), code); true)
    end
end

# The README's table is rows of the one in docs/src/verification.md, which is
# where a re-measurement goes. Each README row must carry that table's tolerance
# and measured value, so the two cannot drift apart; the claim may be worded for
# the README's shorter context.
@testset "The README's verification rows are the manual's" begin
    rows(file) = [strip.(split(chop(strip(l); head = 1, tail = 1), '|'))
                  for l in split(replace(read(file, String), "\r\n" => "\n"), '\n')
                  if startswith(l, "| ") && !startswith(l, "| claim |")]
    root = joinpath(@__DIR__, "..")
    readme = rows(joinpath(root, "README.md"))
    manual = [r[2:end] for r in rows(joinpath(root, "docs", "src", "verification.md"))]
    @test !isempty(readme)
    for r in readme
        @test r[2:end] in manual
    end
end
