# Vasilek — Vlasov Adaptive Simulator of pLasma Electrodynamics and Kinetics

[![CI](https://github.com/korzhimanov/Vasilek.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/korzhimanov/Vasilek.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/korzhimanov/Vasilek.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korzhimanov/Vasilek.jl)

An ongoing project on developing a parallel 2D2P Maxwell — Vlasov — Boltzmann solver on adaptive meshes.

As for now, the following functionality has been implemented:
* PFC scheme on 1D non-uniform grid
* Strang splitting for 1D1V simulations
* 1D Poisson Fourier solver
* 1D FDTD Maxwell solver with PML
* Advection solvers: upwind, Lax—Wendroff, Godunov, semi-Lagrangian
* BGK collision operator

The program has been verified on the following tests:
* [Long-lasting plasma oscillations](verification/plasma-oscillations-1d1v.html)
* [Landau damping](verification/landau-damping-1d1v.html)

Note that these two documents were generated in 2021 and have not been
regenerated since; they are not currently rebuilt by CI.

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
