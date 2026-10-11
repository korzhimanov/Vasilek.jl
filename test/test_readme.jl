using Vasilek

# The README's two usage examples, executed, so that they cannot drift from the
# API. Their number is asserted too, so that a block is not dropped or added
# unnoticed.
@testset "The README examples run" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    blocks = [m.captures[1] for m in eachmatch(r"```julia\r?\n(.*?)```"s, readme)]
    println("  julia blocks found in README.md: ", length(blocks))
    @test length(blocks) == 2
    for code in blocks
        # A fresh module per block, so one cannot lean on names another defined.
        @test (include_string(Module(:READMEExample), code); true)
    end
end
