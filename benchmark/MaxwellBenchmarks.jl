module MaxwellBenchmarks

using BenchmarkTools
using Vasilek: PoissonFourier1D, FDTD1D

const SUITE = BenchmarkGroup()
SUITE["poisson"] = BenchmarkGroup()
SUITE["fdtd"] = BenchmarkGroup()

const SIZES = (100, 1000, 10000)
const Δx = 0.01
const Δt = 0.8*Δx

for N in SIZES
    ρ = [sin(2π*j*Δx) for j = 0:N]
    e = similar(ρ)
    solve! = PoissonFourier1D.generate_solver(ρ, Δx)
    SUITE["poisson"]["solve $N"] = @benchmarkable $solve!($e, $ρ)

    mesh = FDTD1D.YeeMesh1D{Float64}(N)
    pulse = (y = (t, x) -> 0.0, z = (t, x) -> 0.0)
    advance! = FDTD1D.make_advance_fields(mesh, Δt/Δx, pulse, Δt, Δx, 0.0,
                                          FDTD1D.PML(; N = 0, σ_max = 1.0, Δx = Δx, Δt = Δt))
    j = (y = zeros(N + 1), z = zeros(N + 1))
    t = 0.0
    SUITE["fdtd"]["advance $N"] = @benchmarkable $advance!($t, $j)
end

end # module
