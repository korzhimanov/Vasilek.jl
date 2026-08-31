using Vasilek
using Test
using LinearAlgebra

@testset "Test everything" begin
    include("test_aqua.jl")
    include("test_readme.jl")
    include("test_golden.jl")
    include("test_convergence.jl")
    include("test_invariants.jl")
    include("test_symmetry.jl")
    include("test_contracts.jl")
    include("test_threading.jl")
    include("test_allocations.jl")
    include("test_verification.jl")
    include("test_maxwell_solvers.jl")
    include("test_vlasov_solvers.jl")
    include("test_boltzmann_solvers.jl")
end
