module VlasovBenchmarks

using BenchmarkTools

const SUITE = BenchmarkGroup()

include(joinpath(dirname(@__FILE__),"..","src","VlasovSolver","LaxWendroff.jl"))
import .LaxWendroff

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
    advect![N] = LaxWendroff.generate_solver(f₀[N], f[N], v[N]*Δt[N]/Δx[N])

    SUITE["advect $N"] = @benchmarkable advect![$N]()
end

end # module
