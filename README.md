# Vasilek — Vlasov Adaptive Simulator of pLasma Electrodynamics and Kinetics

[![CI](https://github.com/korzhimanov/Vasilek.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/korzhimanov/Vasilek.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/korzhimanov/Vasilek.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korzhimanov/Vasilek.jl)

An ongoing project on developing a parallel 2D2P Maxwell — Vlasov — Boltzmann solver on adaptive meshes.

As for now, the following functionality has been implemented:
* Advection schemes as dispatchable types: upwind, Lax—Wendroff, Godunov
  (piecewise-constant or -linear, with flux limiters), semi-Lagrangian
  (linear, quadratic or cubic B-splines), and PFC on uniform and non-uniform grids
* Strang splitting for 1D1V simulations
* 1D Poisson Fourier solver
* 1D FDTD Maxwell solver with PML
* BGK collision operator

## Verification

Three runnable studies live in `verification/`. They execute directly and write
their figures beside themselves, and they are written in Literate.jl comment
form so they can also be rendered:

```bash
julia --project=verification verification/landau-damping-1d1v.jl
julia --project=verification verification/plasma-oscillations-1d1v.jl
julia --project=verification verification/wakefield.jl
```

Their headline claims are asserted by the test suite rather than left in prose:

```bash
VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

| claim | asserted |
|---|---|
| Landau damping at k = 0.5 | fitted γ within 5% of the tabulated 0.1533 (measured 0.1498, 2.3%) |
| plasma oscillations, uniform grid | \|Δε/ε\| < 0.5% at t = 3000 (measured 0.38%) |
| plasma oscillations, non-uniform grid | \|Δε/ε\| < 6% at t = 3000 (measured 4.85%) |

The wakefield example runs and is stable, but has **no ponderomotive
coupling**: the laser does not drive the wake. See the note in the script.

Unit conventions are in [docs/normalization.md](docs/normalization.md).


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

A scheme value holds no data, so one can be shared across every line of a
multidimensional sweep and across tasks; the workspace is what belongs to the
task. `dest` and `src` must be distinct. Upgrading from 0.1: see
[docs/migration-0.2.md](docs/migration-0.2.md).

This block is executed by the test suite, so it cannot drift from the API.

## Development

```
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Benchmarks live in their own environment:

```
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark benchmark/runbenchmarks.jl
```

See [CHANGELOG.md](CHANGELOG.md) for recent changes, including several
numerically breaking fixes.
