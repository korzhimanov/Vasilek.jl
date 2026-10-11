# Vasilek — Vlasov Adaptive Simulator of pLasma Electrodynamics and Kinetics

[![CI](https://github.com/korzhimanov/Vasilek.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/korzhimanov/Vasilek.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/korzhimanov/Vasilek.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korzhimanov/Vasilek.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://korzhimanov.github.io/Vasilek.jl/dev/)

An ongoing project on developing a parallel 2D2P Maxwell — Vlasov — Boltzmann
solver on adaptive meshes. What exists today is 1D1V on static grids, uniform
in x and uniform or non-uniform in v; 2D2P and adaptive meshes are the goal, not
yet the code.

As for now, the following functionality has been implemented:
* Advection schemes as dispatchable types: upwind, Lax—Wendroff, Godunov
  (piecewise-constant or -linear, with flux limiters), semi-Lagrangian
  (linear, quadratic or cubic B-splines), and PFC on uniform and non-uniform grids
* Strang splitting for 1D1V simulations, and `vlasov_poisson`, the 1D1V
  electrostatic driver (optionally with a collision operator)
* 1D Poisson Fourier solver
* 1D FDTD Maxwell solver with PML
* BGK collision operator

## Verification

The studies in `verification/` hold the solver to analytic theory and to the
literature, and the test suite asserts their claims; CI runs those assertions on
every pull request. The headline rows:

| claim | asserted | measured |
|---|---|---|
| Landau damping rate, k = 0.3, 0.4, 0.5 | γ within 3% of the root of the kinetic dispersion relation | 0.70%, 1.14%, 0.95% |
| Landau real frequency, same three k | ω within 1% of the same root | 0.08%, 0.22%, 0.26% |
| the Vlasov–Poisson flow runs backwards | round-trip error after v → −v falls by at least 3.5× when the grid is halved, at α = 0.05 | 2.82e-3 → 6.48e-4, ×4.3: second order, the limiter clipping the peak |
| strong Landau damping, α = 0.5 | γ₁ within 0.27…0.30, around the literature's 0.281…0.292, which it quotes as −γ; γ₂ within 0.070…0.090, just under the cited 0.0815…0.0858, and refining moves it toward them | 0.2863, 0.0787 (0.0813 at twice the resolution) |
| collisional Landau damping, `BGK` at ν = 0.1, 0.3, 1 | γ within 1.5% and ω within 0.3% of the collisional root on the grid's field, and γ *falling* as ν rises | +0.37%, +0.49%, +0.76%; +0.05%, +0.08%, +0.03%; 0.1547 → 0.1341 → 0.1074 → 0.0660 |
| plasma oscillation frequency | ω within 0.2% of Bohm–Gross √(1+3k²), and the cold ωₚ excluded | 0.018%, against 0.57% for cold |
| two-stream growth rate, kv₀ = 0.4, 0.6, 0.8 | γ within 3% of the warm kinetic root, and the peak in the right place | 0.63%, 1.95%, 0.09% |
| bump-on-tail growth, Arber and Vann's beam | γ and ω within 0.05% and 0.01% of the kinetic root on the grid's field, 0.5% and 0.05% of the continuum's; the wave travels with the beam | 0.016%, 0.003%; 0.131%, 0.014% |
| a plasma echo, field on, is second-order kinetic theory's | pointwise within 1% of `echo_second_order`, the peak within 1% of its size and two steps of its time | 0.48%, 0.48%, t = 29.05 on both |
| a laser wakefield's size is the one linear theory gives | peak within 10% of `linear_wake`, pointwise rms under 15% | 6.4%, 6.4% |

The full table, the nine studies with their commands, and the measurements
behind each claim are in [docs/src/verification.md](docs/src/verification.md),
which the [documentation](https://korzhimanov.github.io/Vasilek.jl/dev/) renders
with each study's output and figures. Unit conventions are in
[docs/src/normalization.md](docs/src/normalization.md).

## Usage

```julia
using Vasilek

src = [1.0 + 0.5*sinpi(2*(i-1)/128) for i = 1:128]
dest = similar(src)
courant = 0.4

scheme = SemiLagrangian(CubicSpline())
ws = workspace(scheme, length(src))     # once, per task; the size must match
advect!(dest, src, scheme, courant, ws)
```

A scheme value holds no buffers, so one can be shared across every line of a
multidimensional sweep and across tasks; the workspace is what belongs to the
task. (`PFCNonUniform` holds its grid's cell widths, so it fits lines of that
grid only.) `dest` and `src` must be distinct, and `|courant| ≤ 1` for every scheme
but `SemiLagrangian`, which has no Courant limit: `advect!` refuses a step past
it rather than return an answer that looks right and is unstable. Upgrading from
0.1: see [docs/src/migration-0.2.md](docs/src/migration-0.2.md).

A whole 1D1V run, linear Landau damping at `k = 0.5`, through the driver the
verification studies use:

```julia
using Vasilek

k = 0.5
x = collect(range(2π/k/32; step = 2π/k/32, length = 32))   # periodic, 32 cells
v = collect(-4.0:0.2:4.0)
f₀ = [exp(-u^2/2)/sqrt(2π)*(1 + 0.01*cos(k*y)) for u in v, y in x]   # f[v, x]
t = collect(0.0:0.1:20.0)

r = vlasov_poisson(x, v, f₀, t; modes = (k,))
r.ε_e        # electric energy per step; r.E_modes[:, 1] the k = 0.5 field
```

The schemes default to `PFCNonUniform` bounded by `f₀`; `scheme_x`, `scheme_v`,
`collisions = BGK(τ)` and the diagnostics are in its docstring.

The public API is the top-level `export` list plus the two named submodules,
`Vasilek.FDTD1D` (`YeeMesh1D`, `PML`, `Yee1D`, `advance!`) and
`Vasilek.PoissonFourier1D` (`PoissonFFT1D`, `solve!`), which are not
re-exported.

These blocks are executed by the test suite, so they cannot drift from the API.

## Development

```
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=verification -e 'using Pkg; Pkg.instantiate()'             # once
julia --project=docs -e 'using Pkg; Pkg.instantiate()'                     # once
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'                # once
julia --project=. -e 'using Pkg; Pkg.test(test_args = ["readme", "docs"])'  # only those files
VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'             # the verification suite: slow, CI runs it on every PR
julia --project=verification verification/two-stream.jl                     # or any other study
julia --project=docs docs/make.jl                                           # VASILEK_DOCS_EXECUTE=1 runs the studies too
julia --project=benchmark benchmark/runbenchmarks.jl                        # advisory
julia --project=benchmark benchmark/workprecision.jl                        # advisory
```

The `verification`, `benchmark` and `docs` environments take the package from
this checkout through `[sources]`, so they need Julia ≥ 1.11 (a docs build that
runs the studies loads theirs from `verification`), and the benchmark baseline
in `benchmark/results.json` was recorded on the development machine, so
`runbenchmarks.jl --strict` means something only there.

See [CHANGELOG.md](CHANGELOG.md) for the changes.
