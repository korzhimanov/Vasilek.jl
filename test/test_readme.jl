using Vasilek

# The README's usage example, executed.
#
# It had drifted: the block referred to `src`, `dest` and `courant` without
# defining any of them, so it could not have run as written, and the file
# carried two near-identical `## Usage` sections. Neither would have been
# noticed by a test suite that never reads the README.
@testset "The README examples run" begin
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    blocks = [m.captures[1] for m in eachmatch(r"```julia\r?\n(.*?)```"s, readme)]
    println("  julia blocks found in README.md: ", length(blocks))
    @test !isempty(blocks)
    for code in blocks
        # A fresh module per block, so one cannot lean on names another defined.
        @test (include_string(Module(:READMEExample), code); true)
    end
end
