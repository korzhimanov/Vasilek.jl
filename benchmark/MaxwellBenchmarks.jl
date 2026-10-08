module MaxwellBenchmarks

using BenchmarkTools
using Vasilek: PoissonFourier1D, FDTD1D, workspace

const SUITE = BenchmarkGroup()
SUITE["poisson"] = BenchmarkGroup()
SUITE["fdtd"] = BenchmarkGroup()

const SIZES = (100, 1000, 10000)
const Δx = 0.01
const Δt = 0.8*Δx

for N in SIZES
    # N points, not N + 1: the FFT lengths were 101 (prime), 1001 and 10001,
    # which timed FFTW's slow path rather than the solver
    ρ = [sin(2π*j*Δx) for j = 0:N-1]
    e = similar(ρ)
    p = PoissonFourier1D.PoissonFFT1D(length(ρ), Δx); ws = workspace(p)
    SUITE["poisson"]["solve $N"] = @benchmarkable PoissonFourier1D.solve!($e, $ρ, $p, $ws)

    mesh = FDTD1D.YeeMesh1D{Float64}(N)
    pulse = (y = (t, x) -> 0.0, z = (t, x) -> 0.0)
    op = FDTD1D.Yee1D(; Δx, Δt, source = pulse,
                      pml = FDTD1D.PML(; N = 0, σ_max = 1.0, Δx = Δx, Δt = Δt))
    j = (y = zeros(N + 1), z = zeros(N + 1))
    t = 0.0
    SUITE["fdtd"]["advance $N"] = @benchmarkable FDTD1D.advance!($mesh, $op, $t, $j)
end

end # module
