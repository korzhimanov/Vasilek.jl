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

Seven runnable studies live in `verification/`. They execute directly and write
their figures beside themselves, and they are written in Literate.jl comment
form so they can also be rendered:

```bash
julia --project=verification verification/landau-damping-1d1v.jl
julia --project=verification verification/plasma-oscillations-1d1v.jl
julia --project=verification verification/wakefield.jl
julia --project=verification verification/two-stream.jl
julia --project=verification verification/bump-on-tail.jl
julia --project=verification verification/plasma-echo.jl
julia --project=verification verification/bgk-equilibrium.jl
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
three constants typed into the test file — the growth rate of *warm*
counter-streaming beams, which is what the runs contain, and the growing root of
the bump-on-tail distribution, on the continuum and on the grid's own
centred-difference field. The cold two-stream closed form remains as its
zero-temperature limit and is checked as such.

| claim | asserted | measured |
|---|---|---|
| Landau damping rate, k = 0.3, 0.4, 0.5 | γ within 3% of the root of the kinetic dispersion relation | 0.71%, 1.14%, 0.94% |
| Landau real frequency, same three k | ω within 1% of the same root | 0.08%, 0.22%, 0.26% |
| and the runs are linear enough for that to be the right target | bounce phase under 2 in every fitting window | 1.70, 0.45, 0.21 |
| the damping fit does not depend on its window | two windows agree within 3% | 0.18%, 0.05%, 0.04% |
| nor does the frequency fit | two windows agree within 1% | 0.066%, 0.000%, 0.300% |
| the Vlasov–Poisson flow runs backwards | round-trip error after v → −v falls by at least 3.5× when the grid is halved, at α = 0.05 | 2.82e-3 → 6.48e-4, ×4.4: second order, the limiter clipping the peak |
| unless the flow has out-run the grid | at α = 0.5 the same refinement recovers nothing | 0.164 → 0.156, ×1.05 |
| and reversibility is bought with positivity | the two most reversible schemes are the two that drive f negative | 0.085 and 0.106 at f = −0.094 and −0.058, against PFC's 0.164 at f ≥ 0 |
| strong Landau damping, α = 0.5 | γ₁ in the literature's −0.281…−0.292; γ₂ within 0.070…0.090, just under the cited 0.0815…0.0858, and refining moves it toward them | 0.2863, 0.0789 (0.0814 at twice the resolution) |
| and the non-uniform velocity grid agrees | both rates within 2% of the uniform grid at the same Δt | 0.05%, 0.6% |
| a drifting plasma damps the same way | boosted mode matches the rest-frame one once the Doppler phase is removed at each sample's own time, and the fitted rates agree within 0.1% | 1.7e-3, 1.8e-3 rad, 0.014% |
| trapping stops the damping on the bounce time | ω_B·t₀ between 6.5 and 8.5 at four amplitudes, and t₀ ∝ α^(−1/2) | 7.09–8.00, slope −0.556 |
| each mode recurs at its own 2π/(kΔv) | within a plasma period, for the seeded mode and the harmonic it generates | 128.6 vs 125.7, 64.3 vs 62.8 |
| a plasma echo, field off, is the closed form's | pointwise within 0.3% of the peak, the peak within two steps of its own, the sign reversing with the kick | 0.12%, t = 15.28 on both |
| out of a mode no moment could see | the echo over 10⁴ times the seeded mode's density at the kick | 2.4·10⁴ |
| and what a scheme keeps of a filament is what it returns | SemiLagrangian < PFC < LaxWendroff < Upwind; upwind returns under 70% | 0.11%, 1.03%, 4.5%, 45%; 55% |
| the loss is truncation | PFC's error falls at least 5× per halving of Δv | 7.4×, 6.9× |
| a plasma echo, field on, is second-order kinetic theory's | pointwise within 1% of `echo_second_order`, the peak within 1% of its size and two steps of its time | 0.48%, 0.47%, t = 29.05 on both |
| which is not the field-off echo | the closed form's peak over 1/0.6 times the run's and 0.5 later | 2.1×, 1.14 |
| and the residual is the run's | it falls 2.5× from Nx = 64 to 128 | 3.4× |
| a nonlinear equilibrium stays put, two thirds of it trapped | f within 0.5% of its peak and the field within 1% through t = 50 | 0.15%, 0.39% |
| and what moves it is the scheme | L² falls and entropy rises; f and field converge at third order, 5× and 4× per halving | 7.6×, 6.1× |
| with a kink in F on the separatrix, the error sits on it | the worst cell within one of the separatrix at two resolutions, converging at first order | on it both times, 1.75× |
| and the equilibrium is this one | ions built on the Poisson sign the docs used to give hold the reversed field; a potential 10% off the ions' drifts 10× more | E/E₀ = −1.000, 16%; 24×, 27× |
| a trapped population above f = 1 holds too | within 3% of its peak, now that PFC's bound comes from f₀ | 0.99%, against 43.6% at the old bound |
| and the old bound is refused rather than run | at fmax = 1 `PFCNonUniform` throws on its first call; unchecked, it takes f over 30% of its peak away | first call; 43.6% |
| plasma oscillation frequency | ω within 0.2% of Bohm–Gross √(1+3k²), and the cold ωₚ excluded | 0.018%, against 0.57% for cold |
| nor does that frequency depend on its window | two windows agree within 0.2% | 0.006% |
| two-stream growth rate, kv₀ = 0.4, 0.6, 0.8 | γ within 3% of the warm kinetic root, and the peak in the right place | 0.63%, 1.95%, 0.09% |
| the same beams at twice the temperature | γ(vt = 0.6) below γ(vt = 0.3), each within 3% of its own warm root | −4.4% measured against −4.4% predicted |
| two-stream at the cold boundary, kv₀ = 1.0 | γ within 8% of the warm root, where the cold form gives exactly zero | 3.12% |
| two-stream stability boundary, kv₀ = 1.2, 1.6 | no growth, where cold and warm theory both give γ = 0 | decays to 0.052, 0.000 |
| bump-on-tail growth, Arber and Vann's beam | γ and ω within 0.05% and 0.01% of the kinetic root on the grid's field, 0.5% and 0.05% of the continuum's; the wave travels with the beam | 0.016%, 0.003%; 0.131%, 0.014% |
| it saturates by trapping, whatever the seed | ω_B/γ between 1.8 and 2.1 at the peak; a seed 1000× larger peaks within 0.5%, ln(1000)/γ earlier | 1.951; 0.04%, 35.00 vs 34.91 |
| and the trapped beam carries the field | the amplitude swings with a period of 1 to 1.6 bounce periods; the flank's slope in ⟨f⟩ falls at least 3× | 1.316; 20× |
| plasma oscillations, uniform grid | \|Δε/ε\| < 0.5% at t = 3000 | 0.38% |
| plasma oscillations, non-uniform grid | \|Δε/ε\| < 6% at t = 3000 | 4.75% |
| and the energy is one instant's | on the uniform grid ε moves by under 3e-4 of itself within a plasma period | 6.0e-5, where summing the kinetic energy after the kick swung it by 1.2e-3 |
| the laser wakefield is a plasma wave | ω within 2% of Bohm–Gross √(ωₚ² + 3Tk²) | 0.23% |
| the driver travels at the grid's group velocity | pulse speed within 2% of `vg_pulse`, at two resolutions | 0.73%, 1.51% |
| and refining moves it toward the continuum | \|v − √(1−n)\| falls when Δx is halved | 0.9487 − 0.8858 → − 0.9261 |
| its wavelength is the driver's | λ within 3% of 2π√(v² − 3T)/ωₚ, for the measured *and* the predicted driver | 0.37%, 1.12% |
| it is phase-locked to the pulse | ω/k within 2% of the pulse velocity, measured independently | 0.58% |
| its size is the one linear theory gives | peak within 10% of `linear_wake`, pointwise rms under 15% | 6.4%, 6.4% |
| and it is the laser that made it | amplitude ∝ a₀² within 8%; ≥10× the unlit control; ≥4× behind the pulse over ahead | 3.3%, 36×, 9.4× |

The echo is the one place a kinetic theory beyond linear order is held to a
number. With the field off it has a closed form, exact in both amplitudes; with
the field on, `echo_second_order` in `test/echo.jl` composes three linear
responses of the Maxwellian — the seed screened, the kick screened, and the
echo's own density polarising the plasma — and the run reproduces the result,
the ringing after the peak included. Both runs are cheap, and both measure
something no moment of `f` does: what a scheme keeps of filaments finer than any
it can show.

The equilibrium study is the converse of every other run: it starts on a
nonlinear stationary state, two thirds of it trapped particles, and asks how
little the solver moves it. The smooth equilibrium's error spreads over phase
space and converges at third order; give the trapped particles a temperature of
their own, which puts a kink in `F` on the separatrix, and the error sits on the
separatrix, twelve times larger, and converges at first order.

The bump-on-tail study is the first run held to numbers past its linear phase.
The growth rate lands on the kinetic root to 0.13%, and on the root of the grid's
own field — the Poisson solve's centred difference returns `sin(kΔx)/(kΔx)` of a
mode's field — to 0.016%. The wave then grows until it traps the beam feeding it,
at a bounce frequency of 1.95 times the growth rate whatever the seed; the
trapped beam swings round the well and the field swings with it; and the
averaged distribution loses the slope the growth ran on.

An eighth study compares the advection schemes on the physics rather than on a
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
task. `dest` and `src` must be distinct, and `|courant| ≤ 1` for every scheme
but `SemiLagrangian`, which has no Courant limit: `advect!` refuses a step past
it rather than return an answer that looks right and is unstable. Upgrading from
0.1: see [docs/migration-0.2.md](docs/migration-0.2.md).

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
