# The canonical scheme list and stepping helper, shared by the golden,
# convergence, invariant and allocation suites so the four cannot drift apart.

using Vasilek

"""
    march!(f, scheme, c, nsteps; ws = workspace(scheme, length(f)))

`nsteps` steps of `scheme` starting from `f`, returning the result. Allocates
one destination buffer and one workspace, then reuses both.
"""
function march!(f₀, scheme, c, nsteps; ws = workspace(scheme, length(f₀)))
    src = copy(f₀)
    dst = similar(src)
    for _ = 1:nsteps
        advect!(dst, src, scheme, c, ws)
        copyto!(src, dst)
    end
    return src
end

"Schemes on a uniform grid, keyed by the names used in test/data/golden.txt."
uniform_schemes(; fmin = 0.0, fmax = 2.0) = [
    ("Upwind",                   Upwind()),
    ("LaxWendroff",              LaxWendroff()),
    ("Godunov_constant",         Godunov(PiecewiseConstant())),
    ("Godunov_linear",           Godunov(PiecewiseLinear())),
    ("Godunov_linear_VanLeer",   Godunov(PiecewiseLinear(), VanLeer())),
    ("SemiLagrangian_linear",    SemiLagrangian(LinearSpline())),
    ("SemiLagrangian_quadratic", SemiLagrangian(QuadraticSpline())),
    ("SemiLagrangian_cubic",     SemiLagrangian(CubicSpline())),
    ("PFC",                      PFC(fmin = fmin, fmax = fmax)),
]
