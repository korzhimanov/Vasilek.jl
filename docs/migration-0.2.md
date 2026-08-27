# Migrating to 0.2

Schemes are values now, and `advect!` dispatches on them. `generate_solver` is
gone, with no compatibility shim: the module names it lived under (`Upwind`,
`PFC`, …) are the new type names, so the old and new API cannot coexist in one
namespace. This was the one point where the plan for this change had to give
way to the language.

## Advection

```julia
# 0.1
advect! = Upwind.generate_solver(f₀, f)          # or (f₀, f, c) to bake c in
advect!(c)

# 0.2
advect!(f, f₀, Upwind(), c)
```

The scheme holds no arrays. One value can be applied to any data, of any size,
from any number of tasks at once — which is the point: a 2D2P step sweeps
O(N²) independent lines per direction, and the old closures captured a single
shared scratch buffer, so they could not be threaded over them.

| 0.1 | 0.2 |
|---|---|
| `Upwind.generate_solver(f₀, f)` | `Upwind()` |
| `LaxWendroff.generate_solver(f₀, f)` | `LaxWendroff()` |
| `Godunov.generate_solver(f₀, f, :Riemann_constant)` | `Godunov(PiecewiseConstant())` |
| `Godunov.generate_solver(f₀, f, :Riemann_linear)` | `Godunov(PiecewiseLinear())` |
| `Godunov.generate_solver(f₀, f, :Riemann_linear; flux_limiter = :VanLeer)` | `Godunov(PiecewiseLinear(), VanLeer())` |
| `SemiLagrangian.generate_solver(f₀, f; interpolation_order = :Cubic)` | `SemiLagrangian(CubicSpline())` |
| `PFC.generate_solver(f₀, f; fₘᵢₙ = a, fₘₐₓ = b)` | `PFC(fmin = a, fmax = b)` |
| `PFCNonUniform.make_advect_1D!(Δx; fₘᵢₙ = a, fₘₐₓ = b)` | `PFCNonUniform(Δx; fmin = a, fmax = b)` |

The `Symbol` options are types, so a typo is a `MethodError` where it is
written rather than a call on `nothing` from inside the hot loop.

## Workspaces

`SemiLagrangian` and `PFCNonUniform` need scratch memory:

```julia
scheme = SemiLagrangian(CubicSpline())
ws = workspace(scheme, length(f))       # once, per task
for _ in 1:nsteps
    advect!(f, f₀, scheme, c, ws)
    copyto!(f₀, f)
end
```

The four-argument `advect!(dest, src, scheme, c)` allocates a workspace on
every call, so pass one explicitly in any loop.

`workspace` returns `nothing` for schemes that need none, so generic code can
always call the five-argument form.

## Note on the fourth argument

For the uniform-grid schemes it is the Courant number `vΔt/Δx`. For
`PFCNonUniform` it is the displacement `vΔt`, a length: a non-uniform grid has
no single Courant number to quote.

## Collisions

```julia
# 0.1
relax! = BGK.generate_solver(f, v, Δt, τ)   # mutated f in place
relax!()

# 0.2
op = BGK(τ)
ws = workspace(op, length(f))
collide!(f, f₀, op, v, Δt, ws)
```

`Landau1P` follows the same shape but is **not exported** and remains
experimental — its closure is inconsistent and the collision integral does not
converge under grid refinement. See its docstring.

## Numerics

Nothing moved. Every scheme is bit-for-bit identical to 0.1, verified in both
directions, and the golden values recorded before this change still match.
