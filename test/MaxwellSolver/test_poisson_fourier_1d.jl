# include(joinpath("..","..","srс","MaxwellSolver","PoissonFourier1D.jl"))
include(joinpath("..","..","src","MaxwellSolver","PoissonFourier1D.jl"))
import .PoissonFourier1D

@test solve_poisson!() == ""
