module Vasilek

include(joinpath("MaxwellSolver", "FDTD1D.jl"))
include(joinpath("MaxwellSolver", "PoissonFourier1D.jl"))

include(joinpath("VlasovSolver", "LaxWendroff.jl"))
include(joinpath("VlasovSolver", "Upwind.jl"))
include(joinpath("VlasovSolver", "Godunov.jl"))
include(joinpath("VlasovSolver", "SemiLagrangian.jl"))
include(joinpath("VlasovSolver", "PFC.jl"))
include(joinpath("VlasovSolver", "PFCNonUniform.jl"))

include(joinpath("BoltzmannSolver", "BGK.jl"))
include(joinpath("BoltzmannSolver", "Landau1P.jl"))

end  # module Vasilek
