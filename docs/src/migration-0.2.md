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

The public API is what `Vasilek` exports, plus the named submodules
`Vasilek.Grids`, `Vasilek.FDTD1D` and `Vasilek.PoissonFourier1D`, whose names
are not re-exported but for the grid types: write `using Vasilek.FDTD1D` for
their names, or qualify them.

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
`[fmin, fmax]`, and throw a `DomainError` naming the cell that does not. In 0.1
`PFC` checked once, at construction, and `PFCNonUniform` never did: data outside
its bounds gave a wrong answer with no error, and now stops at the first call.
A checked call whose `dest` is of a type with no value in `[fmin, fmax]`, a
Float32 `dest` under `fmin = fmax = 0.1`, throws an `ArgumentError` before it
writes `dest`. `vlasov_poisson` refuses a given checked one whose bounds do not
hold `f₀` with an `ArgumentError` before the first step, where its first call
threw a `DomainError`.

A checked call is not only a check. It takes data past a bound by round-off
(`_roundoff`: 1024 ulps of the data's type, Float64 or finer, at the scale of
the cell past the bound and its neighbours; Float32 and Float16 data are taken
only inside the bounds) and clamps what it returns into `[fmin, fmax]`. The
clamp acts only where rounding left a cell outside the bounds, at the ulp level,
or where the data it took were outside by round-off, even at `c = 0`, and moves
the mass by the ulp it cuts, so `checked = true` and `checked = false` can
differ in the last bits of such a cell. Over 10000 driver steps from a plateau
at the bound meeting rough data, the clamp acted on six cells at most, all in
the first 12 steps, and the mass moved no more than with `checked = false`,
whose steps move it by round-off too. `checked = false` removes the check and
the clamp, at compile time.

The positional constructors, `PFC{T, Checked}(fmin, fmax)` and
`PFCNonUniform{T, Checked}(…)`, round the bounds outward into `T`, as the keyword
forms do; they rounded to nearest. Bounds already in `T` keep their bits, but a
Float32 `PFC` built from Float64 bounds now holds each one a bit further out,
so that it takes data at the bounds it was given.

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
`workspace(scheme, n, Float32)`. So is a `PFCNonUniform` workspace's, and an
`OnGrid`'s, whose sub-steps read the line back from it. Both default to their
grid's type, not the data's: Float64 data under a Float32 grid's scheme needs
`workspace(scheme, n, Float64)`, and a narrower one is refused rather than
rounding the data to it. An `OnGrid`'s is checked against `dest` as well as
`src`, since its sub-steps hold `dest`, and on every call: one that takes a
single step never reads the sub-step line and was exact, and is refused all the
same, so that a run whose displacement grows past a cell does not stop halfway.
Code that passed the grid's default `workspace(og, n)` with wider data has to
pass the data's type.

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

`OnGrid(scheme, grid)` puts a scheme on a grid and takes the displacement for
every scheme: it divides by the spacing for one that takes a Courant number,
refuses that on a non-uniform grid, and splits a step wider than a cell into
equal sub-steps that fit. The grid is a `UniformGrid` or a `NonUniformGrid`
(see "Grids" below), a range, or a vector of nodes, which
`Vasilek.Grids.grid` classifies; `(h, n)` is the uniform grid of `n` cells of
`h`. It is a scheme like the others, with a workspace of its own. It replaces
`VlasovPoisson1D1V.line_advector`, a closure added and removed during 0.2;
`substeps` moved with it to `Vasilek.Advection`.

```julia
og = OnGrid(Upwind(), Δx, length(f₀))           # PFCNonUniform carries its grid: OnGrid(p)
ws = workspace(og, length(f₀))
advect!(f, f₀, og, 0.3Δx, ws)                   # the displacement 0.3Δx, c = 0.3
```

## Grids

Grids are types, in `Vasilek.Grids`: a `UniformGrid`, whose cells share one
spacing, and a `NonUniformGrid`, whose cells each have their own width. Whether
a vector of nodes is uniform is decided where its grid is built, by
`grid(z)`, to its rounding and never past 0.5% of a cell; a range is uniform as
it stands, its `step` the spacing. `OnGrid`, `vlasov_poisson`, `setup` and
`make_poisson` take any of them, and `PFCNonUniform` a grid's widths.

The `OnGrid` added during 0.2 took cell widths, which it checked for
uniformity itself, to rounding with a `scale` beside them. It takes the grid
now, and a vector is its nodes:

| during 0.2 | now |
|---|---|
| `OnGrid(s, fill(h, n))` | `OnGrid(s, h, n)` |
| `OnGrid(s, cell_widths(x); scale = maximum(abs, x))` | `OnGrid(s, x)` |
| `OnGrid(p::PFCNonUniform, Δz)` | `OnGrid(p)`, or `OnGrid(p, NonUniformGrid(; widths = Δz))` |
| `og.h`, `og.n` | `og.grid`: `spacing(og.grid)`, `length(og.grid)`; `Vasilek.Advection.substep_bound(og)` |
| `VlasovPoisson1D1V.cell_widths` | `Vasilek.Grids.cell_widths`, still exported where it was |

A vector of widths given as before is refused where it does not increase, as a
uniform grid's widths never do, with a message that says so. Widths that
increase, as a graded grid's may, are taken for nodes: rewrite every such call.
`OnGrid(og, grid)` keeps `og` on a grid with its cells, to the rounding either
carries, and `OnGrid(og, h, n)` is now one of those rather than an error.

Calls that ran and are refused now, or take other bits:

  * a scheme that takes a Courant number on a grid whose spacings differ by
    more than 0.5% of a cell, which its rounding admitted: a Float16 `v` (8
    ulps of its largest node are 16% of a cell on [−4, 4]), a Float16 range
    included, and Float32 grids that reach 10⁵ cells from the origin, as the
    driver notes' "Grids" tabulates; so is an `x` of them, by the field solve,
    whatever the schemes. The default `PFCNonUniform` runs on such a `v` as
    before;
  * every run on a uniform grid moves at round-off: in Float64 near the
    origin 1.9e-15 of the peak at most on the Landau runs of the driver notes'
    "Grids", over 100 to 400 steps; far from it by the bias the mean spacing
    removes, 6.7e-12 at 1e6; and on a Float32 grid by its nodes' rounding,
    1.2e-8 under Float64 data. Its widths are the grid's spacing in
    every cell, in the driver's sums and the collision operators', where they
    were the nodes' cell widths, and a vector's spacing is its mean, where it
    was its first difference; `cell_widths` of a range is its grid's, its step
    in every cell where it is uniform;
  * a collision operator takes a grid's widths through `GridNodes`, which
    `Collide(op, grid, Δt)` hands it; on a vector of nodes `collide!` still
    takes their own cell widths, so a step written out to match the driver
    passes `Vasilek.Collisions.GridNodes(grid(v))`;
  * a `PFCNonUniform` built on the nodes' `cell_widths` beside a uniform grid is
    refused where they are further from the grid's spacing than 4 ulps of the
    coarser type or 1e-12, whichever is more: on the suite's Float32 test `x`,
    8.3e-7 of a cell; far from the origin, 9.2e-10 at 1e6; on enough cells that
    an ulp of the largest node passes 1e-12 of one, 1.6e-12 at 16384 on
    [0, 4π), where 4096 are 1.8e-13 off and taken. Build it on the grid, `PFCNonUniform(grid(x); …)`, as the
    error says; on the suite's common Float64 grids the nodes' widths are
    within 1.4e-14 of it, and are taken;
  * `vlasov_poisson` on ranges takes their nodes collected and their `step` for
    the spacing, where their widths were range arithmetic: a run on ranges is
    the run on the collected nodes to the bit where the step is their mean
    spacing, as on the suite's grids;
  * a Rational `v` runs on its nodes in Float64, where its widths were Rational,
    its spacing the exact mean rounded once;
  * the default scheme on a uniform grid is `PFC`, where it was `PFCNonUniform`
    on the grid's widths: `prob.ox.scheme` and `prob.ov.scheme` are `PFC`s, and
    a run with them moves at round-off, in Float64 3.2e-15 of the peak at most
    over 400 steps on the Landau runs of the driver notes' "The default
    schemes". Given for both directions on Float64 grids,
    `f -> PFCNonUniform(grid(x); fmin = 0, fmax = maximum(f))` and the same on
    `grid(v)` run as the old defaults to the bit, with `fmax = Inf` in place of
    `maximum(f)` under `collisions`, where the defaults have no upper bound;
  * `vlasov_poisson` with a non-uniform `x` throws the field solve's
    `ArgumentError` whatever `scheme_x`, where `OnGrid`'s came first under a
    scheme that takes a Courant number, and before `f₀`'s own `DomainError`.

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
sx  = OnGrid(Upwind(), x)                       # the nodes; or OnGrid(Upwind(), Δx, 32)
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

The workspace is checked as `advect!`'s is: its buffers of the lengths
`workspace(op, length(v))` gives them, a type no narrower than the data's, and
no memory shared with `dest` or `src`. `workspace(BGK(0.5f0), n)` is Float32;
Float64 data needs `workspace(op, n, Float64)`, which the five-argument
`collide!` builds. Some calls that ran are now refused:

  * `BGK` with `dest`, `src` or `v` not one-based, as `Landau1P` already was;
  * a `dest` or `v` sharing any buffer of the workspace, `Landau1P`'s `df`
    included, where only a `dest` over its `F` was refused;
  * a `dest` sharing memory with `v`: `Landau1P` came back 1.6e-4 off with
    `dest` the nodes themselves;
  * a `dest` sharing memory with `src` other than as `src` itself:
    `collide!(f, f, op, v, Δt, ws)` is the step in place, but a `dest` one node
    on from `src` in the same array came back 1.005 off under `Landau1P`.

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
