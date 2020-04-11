using Plots
using Test

@testset "Test everything" begin
    include("test_maxwell_sovers.jl")
    # include(joinpath("..","srс","MaxwellSolver","PoissonFourier1D.jl"))
    # include(joinpath("..","src","MaxwellSolver","PoissonFourier1D.jl"))
end
