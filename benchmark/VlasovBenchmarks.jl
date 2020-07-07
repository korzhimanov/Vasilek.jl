module VlasovBenchmarks

using BenchmarkTools

const SUITE = BenchmarkGroup()

include(joinpath(dirname(@__FILE__),"..","src","VlasovSolver","LaxWendroff.jl"))
import .LaxWendroff
SUITE["LaxWendroff"] = BenchmarkGroup()
SUITE["LaxWendroff c"] = BenchmarkGroup()

include(joinpath(dirname(@__FILE__),"..","src","VlasovSolver","Upwind.jl"))
import .Upwind
SUITE["Upwind"] = BenchmarkGroup()
SUITE["Upwind с"] = BenchmarkGroup()

include(joinpath(dirname(@__FILE__),"..","src","VlasovSolver","Godunov.jl"))
import .Godunov
SUITE["Godunov constant"] = BenchmarkGroup()
SUITE["Godunov linear"] = BenchmarkGroup()
SUITE["Godunov linear VanLeer"] = BenchmarkGroup()

Δx = Dict()
Δt = Dict()
v = Dict()
f₀ = Dict()
f = Dict()
advect! = Dict()

for N in [100, 1000, 10000]
    Δx[N] = 0.01
    Δt[N] = 0.8*Δx[N]
    v[N] = 1
    f₀[N] = [1.0 + 0.01*sin(2π*i*Δx[N]) for i = 0:N]
    f[N] = similar(f₀[N])
    advect![N] = Dict()

    advect![N][:LaxWendroff] = LaxWendroff.generate_solver(f₀[N], f[N], v[N]*Δt[N]/Δx[N])
    SUITE["LaxWendroff"]["advect $N"] = @benchmarkable advect![$N][:LaxWendroff]()

    advect![N][:LaxWendroff_c] = LaxWendroff.generate_solver(f₀[N], f[N])
    SUITE["LaxWendroff c"]["advect $N"] = @benchmarkable advect![$N][:LaxWendroff_c]($(v[N]*Δt[N]/Δx[N]))

    advect![N][:Upwind] = Upwind.generate_solver(f₀[N], f[N], v[N]*Δt[N]/Δx[N])
    SUITE["Upwind"]["advect $N"] = @benchmarkable advect![$N][:Upwind]()

    advect![N][:Upwind_с] = Upwind.generate_solver(f₀[N], f[N])
    SUITE["Upwind с"]["advect $N"] = @benchmarkable advect![$N][:Upwind_с]($(v[N]*Δt[N]/Δx[N]))

    advect![N][:Godunov_constant] = Godunov.generate_solver(f₀[N], f[N], :Riemann_constant)
    SUITE["Godunov constant"]["advect $N"] = @benchmarkable advect![$N][:Godunov_constant]($(v[N]*Δt[N]/Δx[N]))

    advect![N][:Godunov_linear] = Godunov.generate_solver(f₀[N], f[N], :Riemann_linear)
    SUITE["Godunov linear"]["advect $N"] = @benchmarkable advect![$N][:Godunov_linear]($(v[N]*Δt[N]/Δx[N]))

    advect![N][:Godunov_linear_VanLeer] = Godunov.generate_solver(f₀[N], f[N], :Riemann_linear; flux_limiter = :VanLeer)
    SUITE["Godunov linear VanLeer"]["advect $N"] = @benchmarkable advect![$N][:Godunov_linear_VanLeer]($(v[N]*Δt[N]/Δx[N]))
end

end # module
