using Test
include(joinpath("MaxwellSolver", "Maxwell1D.jl"))
using .Maxwell1D

@test advance_fields!() == ""
