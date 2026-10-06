using Vasilek
using Test
using LinearAlgebra

# Every test file, in the order a full run takes them, by its name without
# `test_` and `.jl`. Each file also runs on its own, so arguments select:
#
#     Pkg.test(test_args = ["golden", "contracts"])
#
# runs those two. An unknown name is an error rather than a run of nothing.
const TEST_FILES = ["aqua", "readme", "docs", "dispersion", "golden", "amplification",
                    "convergence", "comparison", "invariants", "symmetry", "contracts",
                    "threading", "allocations", "driver", "harness", "verification",
                    "maxwell_solvers", "em_plasma", "vlasov_solvers", "boltzmann_solvers"]

const SELECTED_TESTS = Set(ARGS)
let unknown = setdiff(SELECTED_TESTS, TEST_FILES)
    isempty(unknown) || error("no test file for ", join(sort!(collect(unknown)), ", "),
                              "; the names are ", join(TEST_FILES, ", "))
end

@testset "Test everything" begin
    for stem in TEST_FILES
        isempty(SELECTED_TESTS) || stem ∈ SELECTED_TESTS || continue
        include("test_$stem.jl")
    end
end
