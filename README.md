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
literature. Each runs on its own and writes its figures beside itself:

```bash
julia --project=verification -e 'using Pkg; Pkg.instantiate()'   # once, Julia ≥ 1.11
julia --project=verification verification/landau-damping-1d1v.jl  # or any other study there
```

Their headline claims are asserted by the test suite rather than left in prose.
CI runs this on every pull request; locally it is behind an environment variable
so that a default `Pkg.test()` stays instant:

```bash
VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

The [documentation](https://korzhimanov.github.io/Vasilek.jl/dev/) renders each
study with its output and figures. Each claim, the tolerance it is held to and
the value measured are tabulated in
[docs/src/verification.md](docs/src/verification.md), with notes on the
studies that needed more than a number.

Unit conventions are in [docs/src/normalization.md](docs/src/normalization.md).


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

These blocks are executed by the test suite, so they cannot drift from the API.

## Development

```
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Benchmarks live in their own environment:

```
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'   # once, Julia ≥ 1.11
julia --project=benchmark benchmark/runbenchmarks.jl
julia --project=benchmark benchmark/workprecision.jl
```

`runbenchmarks.jl` times each kernel against a stored baseline, measured on the
development machine, so a run anywhere else compares machines as much as code;
it prints the baseline's Julia and BenchmarkTools versions beside its own.
`workprecision.jl` pairs error with the cost of reaching it and prints the
efficiency frontier per problem class — the schemes no other scheme beats on
both axes. Both are advisory and exit 0; the accuracy half of the comparison is
gated in `test/test_comparison.jl`, the timing half is not.

The documentation builds with Documenter, in an environment of its own; the
studies it renders load their packages from the verification environment:

```
julia --project=docs -e 'using Pkg; Pkg.instantiate()'           # once, Julia ≥ 1.11
julia --project=verification -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=docs docs/make.jl                                # into docs/build, studies not run
VASILEK_DOCS_EXECUTE=1 julia --project=docs docs/make.jl         # studies run, as CI does
```

See [CHANGELOG.md](CHANGELOG.md) for recent changes, including several
numerically breaking fixes.
