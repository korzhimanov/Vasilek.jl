module Vasilek

"""
    workspace(x, n[, T])

Scratch memory for the scheme, operator or solver `x` at problem size `n`, or
`nothing` when it needs none. One generic for the whole package; each module
adds its methods.
"""
function workspace end

# First: the schemes and the driver dispatch on its grid types.
include("Grids.jl")
include(joinpath("MaxwellSolver", "FDTD1D.jl"))
include(joinpath("MaxwellSolver", "PoissonFourier1D.jl"))

include(joinpath("VlasovSolver", "Advection.jl"))
# Before the splitting, whose collision hook takes an operator.
include(joinpath("BoltzmannSolver", "Collisions.jl"))
include(joinpath("VlasovSolver", "StrangSplitting.jl"))

include(joinpath("VlasovSolver", "VlasovPoisson1D1V.jl"))

using .Grids: AbstractGrid1D, UniformGrid, NonUniformGrid
using .Advection
using .Collisions
using .VlasovPoisson1D1V: vlasov_poisson

export AbstractAdvection1D, advect!, workspace,
       Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC, PFCNonUniform, OnGrid,
       PiecewiseConstant, PiecewiseLinear, NoLimiter, VanLeer, Superbee,
       LinearSpline, QuadraticSpline, CubicSpline,
       AbstractGrid1D, UniformGrid, NonUniformGrid,
       AbstractCollisionOperator, collide!, BGK,
       vlasov_poisson

end  # module Vasilek
