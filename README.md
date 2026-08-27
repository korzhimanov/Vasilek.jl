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

The program has been verified on the following tests:
* [Long-lasting plasma oscillations](verification/plasma-oscillations-1d1v.html)
* [Landau damping](verification/landau-damping-1d1v.html)

Note that these two documents were generated in 2021 and have not been
regenerated since; they are not currently rebuilt by CI.

## Usage

```julia
using Vasilek

scheme = SemiLagrangian(CubicSpline())
ws = workspace(scheme, length(src))    # once, per task
advect!(dest, src, scheme, courant, ws)
```

Upgrading from 0.1: see [docs/migration-0.2.md](docs/migration-0.2.md).

## Usage

```julia
using Vasilek

scheme = SemiLagrangian(CubicSpline())
ws = workspace(scheme, length(src))     # once, per task
advect!(dest, src, scheme, courant, ws)
```

A scheme value holds no data, so one can be shared across every line of a
multidimensional sweep and across tasks. Upgrading from 0.1: see
[docs/migration-0.2.md](docs/migration-0.2.md).

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
