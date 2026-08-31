@testset "Test flux limiters" begin
    include(joinpath("VlasovSolver", "test_flux_limiters.jl"))
end

@testset "Test Vlasov solvers" begin
    include(joinpath("VlasovSolver", "test_1d_advection.jl"))
    include(joinpath("VlasovSolver", "test_nonuniform_advection.jl"))
    include(joinpath("VlasovSolver", "test_strang_splitting.jl"))
end
