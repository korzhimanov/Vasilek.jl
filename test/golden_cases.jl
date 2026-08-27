# Shared by the golden test and its generator, so the two cannot drift apart.

include(joinpath(@__DIR__, "scheme_cases.jl"))

const GOLDEN_N = 64
const GOLDEN_C = 0.4
const GOLDEN_STEPS = 8

golden_initial() = [1.0 + 0.5*sin(2π*(i-1)/GOLDEN_N) for i = 1:GOLDEN_N]

golden_run(scheme) = march!(golden_initial(), scheme, GOLDEN_C, GOLDEN_STEPS)

const GOLDEN_CASES = uniform_schemes()
