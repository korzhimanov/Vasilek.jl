module VlasovBenchmarks

using BenchmarkTools
using Vasilek

const SUITE = BenchmarkGroup()

const SIZES = (100, 1000, 10000)
const Δx = 0.01
const COURANT = 0.8

initial(N) = [1.0 + 0.01*sin(2π*i*Δx) for i = 0:N]

const SCHEMES = (
    ("LaxWendroff",              LaxWendroff()),
    ("Upwind",                   Upwind()),
    ("Godunov constant",         Godunov(PiecewiseConstant())),
    ("Godunov linear",           Godunov(PiecewiseLinear())),
    ("Godunov linear VanLeer",   Godunov(PiecewiseLinear(), VanLeer())),
    ("SemiLagrangian linear",    SemiLagrangian(LinearSpline())),
    ("SemiLagrangian quadratic", SemiLagrangian(QuadraticSpline())),
    ("SemiLagrangian cubic",     SemiLagrangian(CubicSpline())),
    ("PFC",                      PFC(fmin = 0.0, fmax = 2.0)),
)

for (name, _) in SCHEMES
    SUITE[name] = BenchmarkGroup()
end
SUITE["PFCNonUniform"] = BenchmarkGroup()

for N in SIZES
    src = initial(N)
    n = length(src)

    for (name, scheme) in SCHEMES
        dst = similar(src)
        ws = workspace(scheme, n)
        # Values are interpolated, not looked up: leaving a lookup in a
        # non-const global inside the timed region measured the lookup.
        SUITE[name]["advect $N"] = @benchmarkable advect!($dst, $src, $scheme, $COURANT, $ws)
    end

    nu = PFCNonUniform(fill(Δx, n); fmin = 0.0, fmax = 2.0)
    dst = similar(src)
    ws = workspace(nu, n)
    α = COURANT*Δx
    SUITE["PFCNonUniform"]["advect $N"] = @benchmarkable advect!($dst, $src, $nu, $α, $ws)
end

end # module
