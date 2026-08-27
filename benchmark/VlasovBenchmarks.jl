module VlasovBenchmarks

using BenchmarkTools
using Vasilek: LaxWendroff, Upwind, Godunov, SemiLagrangian, PFC, PFCNonUniform

const SUITE = BenchmarkGroup()

const SIZES = (100, 1000, 10000)
const Δx = 0.01
const COURANT = 0.8

initial(N) = [1.0 + 0.01*sin(2π*i*Δx) for i = 0:N]

# Schemes taking the Courant number at call time -- the path that survives the
# planned move to scheme structs.
const SCHEMES = (
    ("LaxWendroff",            (a, b) -> LaxWendroff.generate_solver(a, b)),
    ("Upwind",                 (a, b) -> Upwind.generate_solver(a, b)),
    ("Godunov constant",       (a, b) -> Godunov.generate_solver(a, b, :Riemann_constant)),
    ("Godunov linear",         (a, b) -> Godunov.generate_solver(a, b, :Riemann_linear)),
    ("Godunov linear VanLeer", (a, b) -> Godunov.generate_solver(a, b, :Riemann_linear; flux_limiter = :VanLeer)),
    ("SemiLagrangian linear",    (a, b) -> SemiLagrangian.generate_solver(a, b; interpolation_order = :Linear)),
    ("SemiLagrangian quadratic", (a, b) -> SemiLagrangian.generate_solver(a, b; interpolation_order = :Quadratic)),
    ("SemiLagrangian cubic",     (a, b) -> SemiLagrangian.generate_solver(a, b; interpolation_order = :Cubic)),
    ("PFC",                    (a, b) -> PFC.generate_solver(a, b; fₘᵢₙ = 0.0, fₘₐₓ = 2.0)),
)

# The same kernels reached through the overload that bakes the Courant number
# in at construction. Kept so the delegation introduced when the two were
# unified can be shown to cost nothing.
const BAKED = (
    ("LaxWendroff baked", (a, b, c) -> LaxWendroff.generate_solver(a, b, c)),
    ("Upwind baked",      (a, b, c) -> Upwind.generate_solver(a, b, c)),
    ("PFC baked",         (a, b, c) -> PFC.generate_solver(a, b, c; fₘᵢₙ = 0.0, fₘₐₓ = 2.0)),
)

for (name, _) in SCHEMES
    SUITE[name] = BenchmarkGroup()
end
for (name, _) in BAKED
    SUITE[name] = BenchmarkGroup()
end
SUITE["PFCNonUniform"] = BenchmarkGroup()

for N in SIZES
    f₀ = initial(N)

    for (name, mk) in SCHEMES
        step! = mk(f₀, similar(f₀))
        # The closure is interpolated, not looked up. Writing `advect![$N][:X]()`
        # left a non-const global Dict lookup inside the timed region, which at
        # N = 100 plausibly cost more than the kernel.
        SUITE[name]["advect $N"] = @benchmarkable $step!($COURANT)
    end

    for (name, mk) in BAKED
        step! = mk(f₀, similar(f₀), COURANT)
        SUITE[name]["advect $N"] = @benchmarkable $step!()
    end

    Δ = fill(Δx, N + 1)
    g = initial(N)
    nustep! = PFCNonUniform.make_advect_1D!(Δ; fₘᵢₙ = 0.0, fₘₐₓ = 2.0)
    α = COURANT*Δx
    SUITE["PFCNonUniform"]["advect $N"] = @benchmarkable $nustep!($g, $α)
end

end # module
