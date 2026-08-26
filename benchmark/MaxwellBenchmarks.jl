module MaxwellBenchmarks

using BenchmarkTools
using Vasilek: PoissonFourier1D

const SUITE = BenchmarkGroup()


Δx = Dict()
ρ = Dict()
e = Dict()
solve! = Dict()

for N in [100, 1000, 10000]
    Δx[N] = 0.01
    ρ[N] = [sin(2π*j*Δx[N]) for j = 0:N]
    e[N] = similar(ρ[N])
    solve![N] = PoissonFourier1D.generate_solver(ρ[N], Δx[N])
    SUITE["poisson $N"] = @benchmarkable solve![$N]($(e[N]), $(ρ[N]))
end

end # module
