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

Four runnable studies live in `verification/`. They execute directly and write
their figures beside themselves, and they are written in Literate.jl comment
form so they can also be rendered:

```bash
julia --project=verification verification/landau-damping-1d1v.jl
julia --project=verification verification/plasma-oscillations-1d1v.jl
julia --project=verification verification/wakefield.jl
julia --project=verification verification/two-stream.jl
```

Their headline claims are asserted by the test suite rather than left in prose.
CI runs this on every pull request; locally it is behind an environment variable
so that a default `Pkg.test()` stays instant:

```bash
VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

The analytic side of those claims is computed rather than quoted:
`test/dispersion.jl` solves the kinetic dispersion relation through the plasma
dispersion function, which gives the Landau roots at any `k` — they used to be
three constants typed into the test file — and the growth rate of *warm*
counter-streaming beams, which is what the runs contain. The cold two-stream
closed form remains as its zero-temperature limit and is checked as such.

| claim | asserted | measured |
|---|---|---|
| Landau damping rate, k = 0.3, 0.4, 0.5 | γ within 3% of the root of the kinetic dispersion relation | 0.71%, 1.09%, 1.11% |
| Landau real frequency, same three k | ω within 1% of the same root | 0.08%, 0.07%, 0.14% |
| and the runs are linear enough for that to be the right target | bounce phase under 2 in every fitting window | 1.70, 0.45, 0.21 |
| the damping fit does not depend on its window | two windows agree within 3% | 0.19%, 0.05%, 0.36% |
| nor does the frequency fit | two windows agree within 1% | 0.066%, 0.146%, 0.100% |
| the Vlasov–Poisson flow runs backwards | round-trip error after v → −v falls by at least 4× when the grid is halved, at α = 0.05 | 1.97e-3 → 3.44e-4, ×5.7 |
| unless the flow has out-run the grid | at α = 0.5 the same refinement recovers nothing | 0.164 → 0.156, ×1.05 |
| and reversibility is bought with positivity | the two most reversible schemes are the two that drive f negative | 0.085 and 0.106 at f = −0.094 and −0.058, against PFC's 0.164 at f ≥ 0 |
| strong Landau damping, α = 0.5 | γ₁ in the literature's −0.281…−0.292 and γ₂ in its 0.0770…0.0858; refining moves γ₂ toward them | 0.2863, 0.0789 |
| and the non-uniform velocity grid agrees | both rates within 2% of the uniform grid at the same Δt | 0.04%, 0.6% |
| a drifting plasma damps the same way | boosted mode matches the rest-frame one once the Doppler phase is removed, and the fitted rates agree | 1.5e-3, 8.0e-3 rad, 0.133% |
| trapping stops the damping on the bounce time | ω_B·t₀ between 6.5 and 8.5 at four amplitudes, and t₀ ∝ α^(−1/2) | 7.09–8.01, slope −0.557 |
| each mode recurs at its own 2π/(kΔv) | within a plasma period, for the seeded mode and the harmonic it generates | 128.6 vs 125.7, 64.3 vs 62.8 |
| plasma oscillation frequency | ω within 0.2% of Bohm–Gross √(1+3k²), and the cold ωₚ excluded | 0.018%, against 0.57% for cold |
| nor does that frequency depend on its window | two windows agree within 0.2% | 0.006% |
| two-stream growth rate, kv₀ = 0.4, 0.6, 0.8 | γ within 3% of the warm kinetic root, and the peak in the right place | 0.39%, 1.95%, 0.10% |
| the same beams at twice the temperature | γ(vt = 0.6) below γ(vt = 0.3), each within 3% of its own warm root | −4.4% measured against −4.4% predicted |
| two-stream at the cold boundary, kv₀ = 1.0 | γ within 8% of the warm root, where the cold form gives exactly zero | 3.19% |
| two-stream stability boundary, kv₀ = 1.2, 1.6 | no growth, where cold and warm theory both give γ = 0 | decays to 0.053, 0.000 |
| plasma oscillations, uniform grid | \|Δε/ε\| < 0.5% at t = 3000 | 0.38% |
| plasma oscillations, non-uniform grid | \|Δε/ε\| < 6% at t = 3000 | 4.85% |
| the laser wakefield is a plasma wave | ω within 2% of Bohm–Gross √(ωₚ² + 3Tk²) | 0.23% |
| the driver travels at the grid's group velocity | pulse speed within 2% of `vg_pulse`, at two resolutions | 0.73%, 1.51% |
| and refining moves it toward the continuum | \|v − √(1−n)\| falls when Δx is halved | 0.9487 − 0.8858 → − 0.9261 |
| its wavelength is the driver's | λ within 3% of 2π√(v² − 3T)/ωₚ, for the measured *and* the predicted driver | 0.37%, 1.12% |
| it is phase-locked to the pulse | ω/k within 2% of the pulse velocity, measured independently | 0.58% |
| its size is the one linear theory gives | peak within 10% of `linear_wake`, pointwise rms under 15% | 6.4%, 6.4% |
| and it is the laser that made it | amplitude ∝ a₀² within 8%; ≥10× the unlit control; ≥4× behind the pulse over ahead | 3.3%, 36×, 9.4× |

A fifth study compares the advection schemes on the physics rather than on a
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
rows above measure the wake it produces against linear wakefield theory rather
than against a bound.

The study runs at twenty cells per laser wavelength, and that is a physics
choice rather than a taste. At the ten it used before, the Yee dispersion
relation puts the driver's group velocity 4.5% below `√(1−n)`: the pulse arrives
late and the wake keeps station with it, so the wake's phase velocity is wrong
by the same amount — `γ_φ ≈ 2.2` where the physics gives 3.0. None of the
assertions about the wake's *frequency* could see that, and the comparison
against `linear_wake` cannot either, since that reference is driven by the `Φ`
of the run it is checking. It took a closed form for the driver, `vg_pulse`, and
a second resolution.

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
