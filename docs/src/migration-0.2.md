# Migrating to 0.2

Schemes are values now, and `advect!` dispatches on them. The 0.1 closure
factories are gone, with no compatibility shim: the advection and collision
`generate_solver`s, the Poisson one, and `FDTD1D.make_advance_fields`. For the
advection schemes there was no choice: the module names `generate_solver` lived
under (`Upwind`, `PFC`, …) are the new type names, so the old and new API cannot
coexist in one namespace. This was the one point where the plan for this change
had to give way to the language. The others follow them so that 0.2 has one
shape throughout, and since no release carries them, a shim would have
preserved nothing.

The public API is what `Vasilek` exports, plus the two named submodules
`Vasilek.FDTD1D` and `Vasilek.PoissonFourier1D`, which are not re-exported:
write `using Vasilek.FDTD1D` for their names, or qualify them.

## Advection

```julia
# 0.1
advect! = Upwind.generate_solver(f₀, f)          # or (f₀, f, c) to bake c in
advect!(c)

# 0.2
advect!(f, f₀, Upwind(), c)
```

The scheme holds no arrays, `PFCNonUniform` aside: it carries its grid's cell
widths and limiter coefficients, so it fits lines of that grid only, though it
is still read-only and shareable between tasks. Every other value can be
applied to any data, of any size, from any number of tasks at once — which is
the point: a 2D2P step sweeps
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

The two `:Riemann_linear` rows change the numbers as well as the name. In 0.1
the flux was the reconstruction's value at the interface, so the update was
forward Euler on a limited slope. That is first order, total-variation
diminishing only to `|c| ≤ 1/2`, and unstable at every `c` without a limiter.
`Godunov(PiecewiseLinear())` averages the reconstruction over the strip that
crosses the interface in one step, which puts a `(1 − |c|)` factor on the slope
term. The result is Sweby's flux-limited Lax–Wendroff: second order,
total-variation diminishing to `|c| = 1` with `VanLeer`, and without a limiter
exactly `LaxWendroff()`. A run that used either row will not reproduce its 0.1
numbers. With `VanLeer` at `c = 0.4` and N = 512, the L² error after one
traversal falls from 7.03e-3 to 8.55e-5 on a sine, and rises from 8.45e-3 to
2.77e-2 on a square pulse, because 0.1's flux also steepened the jump.

`PFC` and `PFCNonUniform` check on every call that the data lie within
`[fmin, fmax]`, and throw a `DomainError` when they do not. In 0.1 `PFC`
checked once, at construction, and `PFCNonUniform` never did: data outside its
bounds gave a wrong answer with no error, and now stops at the first call.
`checked = false` removes the check, at compile time.

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

The size has to match exactly: `workspace(scheme, n)` is what a call on `n`
elements requires, and both a shorter and a longer buffer are rejected. A
longer one is not merely wasteful — the spline prefilter runs over the whole
buffer, so a `SemiLagrangian` handed a workspace built for a longer line used
to return garbage silently. Allocate one per line length, not one big one for
all of them. A `SemiLagrangian` workspace also belongs to its spline: it holds
that spline's prefilter, factorised when it is built, and one built for another
spline is refused. Its element type is the type the step computes in, whatever
the type of `c`: it may be wider than the data, not narrower, and it has to hold
every knot exactly, so `Float16` data on more than 2047 cells needs
`workspace(scheme, n, Float32)`.

`advect!` also requires that `dest` and `src` share no memory (a view of `src`
counts) and at least three cells. Both used to be quietly wrong rather than an
error.

So is a step past a scheme's Courant limit, and it is refused now too: `|c| > 1`
for every scheme but `SemiLagrangian`, and for `PFCNonUniform` a displacement
wider than its narrowest cell, raise a `DomainError`. 0.1 took the step and
returned an answer that looked right for tens of steps before it grew without
bound. A run that needs the longer step can split it into sub-steps that fit,
as `OnGrid` (below) does, or use `SemiLagrangian`.

## Note on the fourth argument

For the uniform-grid schemes it is the Courant number `vΔt/Δx`. For
`PFCNonUniform` it is the displacement `vΔt`, a length: a non-uniform grid has
no single Courant number to quote.

`OnGrid(scheme, Δz)` puts a scheme on the grid of cell widths `Δz` and takes
the displacement for every scheme: it divides by the spacing for a uniform-grid
one, refuses that on a non-uniform grid, and splits a step wider than a cell
into equal sub-steps that fit. It is a scheme like the others, with a workspace
of its own. It replaces `VlasovPoisson1D1V.line_advector`, a closure added and
removed during 0.2; `substeps` moved with it to `Vasilek.Advection`.

```julia
og = OnGrid(Upwind(), fill(Δx, length(f₀)))     # PFCNonUniform carries its grid: OnGrid(p)
ws = workspace(og, length(f₀))
advect!(f, f₀, og, 0.3Δx, ws)                   # the displacement 0.3Δx, c = 0.3
```

## Strang splitting

```julia
# 0.1
g = f'                                          # f[v, x]; g its x lines
make_time_step_2d!((g, f), (_ -> v*Δt, ff -> e*Δt), (advect_x!, advect_v!))

# 0.2
using Vasilek.StrangSplitting: strang_step!, Collide
using Vasilek.VlasovPoisson1D1V: cell_widths
x   = collect(Δx .* (1:32))
fxv = [f₀[j]*(1 + 0.1cos(2π*i/32)) for i in eachindex(x), j in eachindex(v)]   # f[x, v]
sx  = OnGrid(Upwind(), fill(Δx, length(x)))
sv  = OnGrid(PFCNonUniform(cell_widths(v); fmin = 0.0, fmax = Inf))
ws  = workspace(sx, sv, fxv, BGK(τ))            # one per task; the operator's scratch too
E   = 0.1 .* sin.(2π .* (1:32) ./ 32)
strang_step!(fxv, sx, sv, v .* Δt, ft -> E .* Δt, ws, Collide(BGK(τ), v, Δt))
```

`strang_step!(f, sx, sv, cx, cv, ws[, hook])` steps one matrix, `f[x, v]`: the
x sweep runs over its columns and the v sweep over the transposed copy the
workspace keeps, so every line is contiguous. `make_time_step_2d!` is gone,
with its second array and its per-line closures.

  * The schemes are values, on their grids through `OnGrid` when `cx` and `cv`
    are displacements; `cx[j]` is column `j`'s full step, and each half step
    takes `cx[j]/2`.
  * **`cv` is handed the state transposed**, `ft[v, x]`, after the first half
    step, not `f`: a `cv` that summed the rows of `f` sums the columns of `ft`.
    That is the layout the density is summed over, and where the field is
    solved.
  * A collision operator goes in as the hook `Collide(op, v, Δt)`, collided for
    `Δt/2` either side of the kick, `X(Δt/2) · C(Δt/2) K(Δt) C(Δt/2) · X(Δt/2)`,
    with its scratch from `workspace(sx, sv, f, op)`.

## Poisson

```julia
# 0.1
solve! = PoissonFourier1D.generate_solver(ρ₀, Δx)   # removed, with no shim
solve!(e, ρ)

# 0.2
using Vasilek.PoissonFourier1D: PoissonFFT1D, solve!
p = PoissonFFT1D(length(ρ), Δx)          # derivative = :spectral for the exact one
ws = workspace(p)                        # one per task
solve!(e, ρ, p, ws)
```

## FDTD

```julia
# 0.1
mesh = FDTD1D.YeeMesh1D{Float64}(200)
advance_fields! = FDTD1D.make_advance_fields(mesh, Δt/Δx, pulse, Δt, Δx, x_min,
                                             FDTD1D.PML(10, 1e3, Δx, Δt))
advance_fields!(t, j)

# 0.2
using Vasilek.FDTD1D
mesh = YeeMesh1D{Float64}(200)
pulse = (y = (t, x) -> exp(-(x - t + 3)^2), z = (t, x) -> 0.0)   # (t, x) -> amplitude
op = Yee1D(; Δx, Δt, source = pulse)    # cfl = Δt/Δx, x_min = 0.0, pml = PML(; N = 10, σ_max = 1e3, Δx, Δt)
j = (y = zeros(201), z = zeros(201))    # −J·Δt on the N + 1 electric nodes
for s in 1:nsteps
    advance!(mesh, op, s*Δt, j)
end
```

| 0.1 | 0.2 |
|---|---|
| `make_advance_fields(f, cfl, pulse, Δt, Δx, x_min, pml)` | `op = Yee1D(; Δx, Δt, cfl, source = pulse, x_min, pml)` |
| `advance_fields!(t, j)` | `advance!(f, op, t, j)`, which returns `f` |
| `YeeMesh1D{T, S}`, `S` the type of `N` | `YeeMesh1D{T}`, `N` an `Int` |
| `PML{I, T}`, `I` the type of `N`; a `MethodError` for an integer `σ_max` | `PML{T}`; `PML(10, 1000, Δx, Δt)` is `PML(10, 1e3, Δx, Δt)` |

`Yee1D` holds the step's constants, the absorbing layers and the source, and no
fields, so one value steps any number of meshes from any number of tasks, and
`workspace(op)` is `nothing`. The closure captured its mesh; the operator does
not, so the check that the mesh holds both layers and an interior moved from
construction to every `advance!` call, beside a check that `j.y` and `j.z`
have the axes of `mesh.ey`, `1:N + 1`. `cfl` still has to equal `Δt/Δx`, and is
checked when the operator is built, as are steps that are finite and positive
and a `pml` built for this `Δx` and `Δt`, which a `PML` now records. It is a keyword of its own because a caller
that sets the Courant number wants exactly that number in the interior, and
`cfl*Δx/Δx` need not round back to it.

The step works in the type of `Δx` and `Δt`: `cfl` and the layer are converted
to it, so `Float32` steps with a literal `cfl = 0.8` or the default layer give a
`Float32` operator, where the closure computed those terms in `Float64`.
`x_min` keeps its own type.

`FDTD1D` exports `PML` now, beside `YeeMesh1D`, `Yee1D` and `advance!`. The
update is the closure's, moved: with `Float64` data every field comes out
bit-for-bit as `make_advance_fields` left it; with `Float32` data and a
`Float64` layer, to the `Float32` rounding of the fields.

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

`BGK(τ)` now relaxes to the *discrete* Maxwellian, which conserves density,
momentum and energy to round-off on the grid it is given; 0.1's sampled
Maxwellian, which a narrow velocity window cooled, is `BGK(τ; conservative =
false)`.

`Landau1P` follows the same shape but is **not exported**: reach it as
`Vasilek.Collisions.Landau1P`. It is a one-dimensional model of the Landau
operator whose transverse velocities are a bath at `Tₜ`. Its kernel
`2Tₜ/(u² + 2Tₜ)^(3/2)` is finite where 0.1's `2Tₜ/|u|³` was not, so the
collision integral converges under grid refinement; the bath's drag makes the
Maxwellian at `Tₜ` its equilibrium; and the update is in flux form, so the
density is conserved to round-off. The step is forward Euler on a diffusion:
its docstring gives the time-step bound. `Tₜ` defaults to 1, the plasma's own
temperature in thermal units, where it was 1e-3: it no longer scales the rate,
and a bath that cold collapses a thermal line below any practical grid. It
must be positive, and at least the square of the grid's coarsest spacing.

## Numerics

The refactor itself moved nothing: every scheme came out bit-for-bit identical
to 0.1. Several fixes made since do move numbers, on purpose:

  * `Godunov(PiecewiseLinear(), …)`, the `:Riemann_linear` rows above, carries
    the `(1 − |c|)` factor and is a different, second-order scheme;
  * `BGK` relaxes to the discrete Maxwellian by default;
  * the z-polarised source of `FDTD1D` launches its pulse the other way;
  * the quadratic and cubic `SemiLagrangian` solve their periodic prefilter
    themselves rather than through Interpolations.jl, and move in the last bits
    only: 1.2e-15 relative at most over the golden run, under 5e-15 of the
    data's maximum in a wider comparison. The linear spline is unchanged;
  * `Landau1P` has a regularised kernel, a drag, and a flux-form update.

The CHANGELOG lists each with what it changes.
