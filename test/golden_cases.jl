# Shared by the golden test and its generator, so the two cannot drift apart.

const GOLDEN_N = 64
const GOLDEN_C = 0.4
const GOLDEN_STEPS = 8

golden_initial() = [1.0 + 0.5*sin(2π*(i-1)/GOLDEN_N) for i = 1:GOLDEN_N]

"""
    golden_run(mk)

Eight steps of the scheme built by `mk(src, dst)`, feeding each result back
through the captured input buffer.
"""
function golden_run(mk)
    src = golden_initial()
    dst = similar(src)
    step! = mk(src, dst)
    for _ = 1:GOLDEN_STEPS
        step!(GOLDEN_C)
        copyto!(src, dst)
    end
    return src
end

const GOLDEN_CASES = [
    ("Upwind",                   (a,b) -> Upwind.generate_solver(a,b)),
    ("LaxWendroff",              (a,b) -> LaxWendroff.generate_solver(a,b)),
    ("Godunov_constant",         (a,b) -> Godunov.generate_solver(a,b,:Riemann_constant)),
    ("Godunov_linear",           (a,b) -> Godunov.generate_solver(a,b,:Riemann_linear)),
    ("Godunov_linear_VanLeer",   (a,b) -> Godunov.generate_solver(a,b,:Riemann_linear; flux_limiter=:VanLeer)),
    ("SemiLagrangian_linear",    (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Linear)),
    ("SemiLagrangian_quadratic", (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Quadratic)),
    ("SemiLagrangian_cubic",     (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Cubic)),
    ("PFC",                      (a,b) -> PFC.generate_solver(a,b;fₘᵢₙ=0.0,fₘₐₓ=2.0)),
]
