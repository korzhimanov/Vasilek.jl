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

Their headline claims are asserted by the test suite rather than left in prose.
CI runs this on every pull request; locally it is behind an environment variable
so that a default `Pkg.test()` stays instant:

```bash
VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

| claim | asserted | measured |
|---|---|---|
| Landau damping rate, k = 0.3, 0.4, 0.5 | γ within 3% of the tabulated root | 0.42%, 0.51%, 1.45% |
| Landau real frequency, same three k | ω within 1% of the tabulated root | 0.25%, 0.36%, 0.14% |
| the damping fit does not depend on its window | two windows agree within 3% | 1.45%, 0.44%, 0.43% |
| nor does the frequency fit | two windows agree within 1% | 0.05%, 0.15%, 0.10% |
| plasma oscillation frequency | ω within 0.2% of Bohm–Gross √(1+3k²), and the cold ωₚ excluded | 0.018%, against 0.57% for cold |
| nor does that frequency depend on its window | two windows agree within 0.2% | 0.006% |
| two-stream growth rate, kv₀ = 0.4, 0.6, 0.8 | γ within 6% of the closed-form cold root, and the peak in the right place | 1.87%, 3.14%, 0.31% |
| two-stream stability boundary, kv₀ = 1.2, 1.6 | no growth where the closed form gives γ = 0 exactly | decays to 0.053, 0.000 |
| plasma oscillations, uniform grid | \|Δε/ε\| < 0.5% at t = 3000 | 0.38% |
| plasma oscillations, non-uniform grid | \|Δε/ε\| < 6% at t = 3000 | 4.85% |
| the laser wakefield is a plasma wave | ω within 2% of Bohm–Gross √(ωₚ² + 3Tk²) | 0.36% |
| its wavelength is the driver's | λ within 3% of 2π√(v² − 3T)/ωₚ, at the measured pulse speed | 1.16% |
| it is phase-locked to the pulse | ω/k within 2% of the pulse velocity, measured independently | 0.75% |
| its size is the one linear theory gives | peak within 10% of `linear_wake`, pointwise rms under 15% | 4.4%, 8.9% |
| and it is the laser that made it | amplitude ∝ a₀² within 8%; ≥10× the unlit control; ≥4× behind the pulse over ahead | 0.45%, 32×, 7.5× |

A fourth study compares the advection schemes on the physics rather than on a
shifted sine, and is advisory rather than asserted:

```bash
julia --project=verification verification/scheme-comparison.jl
```

Ranked by error in the Landau damping rate, `LaxWendroff` (0.64%) and cubic
`SemiLagrangian` (1.07%) lead and upwind trails at 48.8% — its own dissipation
being two orders of magnitude larger than the damping it is measuring. At 50%
amplitude the two schemes that lead are exactly the two that drive `f` negative,
which is why the solvers default to `PFC`.

The wakefield example was asserted only to run and stay bounded until recently,
because that was all it could support: it had **no ponderomotive coupling**, so
the laser never entered the longitudinal push and the wake it drew was the slab
edges relaxing. Its numbers came out bit-identical whether the transverse
current was right or wrong by thirty-two orders of magnitude. The coupling —
the force `−∂(pʸ² + pᶻ²)/2∂x` in the momentum advection — is now there, and the
five rows above measure the wake it produces against linear wakefield theory
rather than against a bound.

Two approximations remain, and the defaults are chosen to stay inside them
rather than to be impressive: the ponderomotive potential is the
non-relativistic one and the transverse current is taken through momentum
rather than velocity, both of which are corrections of relative order `p⊥²`.
At `a₀ = 0.3` that is 3.2%. The transverse momentum itself is exact — it is the
canonical `p⊥ = −A⊥`, not a force integral. See the docstring on `wakefield`.

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
julia --project=benchmark benchmark/workprecision.jl
```

`runbenchmarks.jl` times each kernel against a stored baseline.
`workprecision.jl` pairs error with the cost of reaching it and prints the
efficiency frontier per problem class — the schemes no other scheme beats on
both axes. Both are advisory and exit 0; the accuracy half of the comparison is
gated in `test/test_comparison.jl`, the timing half is not.

See [CHANGELOG.md](CHANGELOG.md) for recent changes, including several
numerically breaking fixes.
