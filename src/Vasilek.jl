module Vasilek

include(joinpath("MaxwellSolver", "FDTD1D.jl"))
include(joinpath("MaxwellSolver", "PoissonFourier1D.jl"))

include(joinpath("VlasovSolver", "Advection.jl"))
include(joinpath("VlasovSolver", "StrangSplitting.jl"))

include(joinpath("BoltzmannSolver", "Collisions.jl"))

using .Advection
using .Collisions

export AbstractAdvection1D, advect!, workspace,
       Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC, PFCNonUniform,
       PiecewiseConstant, PiecewiseLinear, NoLimiter, VanLeer,
       LinearSpline, QuadraticSpline, CubicSpline,
       AbstractCollisionOperator, collide!, BGK

end  # module Vasilek
