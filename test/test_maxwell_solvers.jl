@testset "Test Maxwell solvers" begin
    include(joinpath("MaxwellSolver", "test_poisson_fourier_1d.jl"))
    include(joinpath("MaxwellSolver", "test_fdtd_1d.jl"))
end
