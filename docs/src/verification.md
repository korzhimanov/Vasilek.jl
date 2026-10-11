# Verification

The solver is held to analytic theory and to the literature by the studies in
[`verification/`](https://github.com/korzhimanov/Vasilek.jl/blob/master/verification),
one script each. Every study has a page of its own under Verification in the
navigation: the script's code, the numbers it prints and the figures it draws.
This page says what the test suite asserts about them, and gathers the
measurements behind each claim, study by study.

## Running the studies

Each study executes directly and writes its figures beside itself. The scripts
are written in Literate.jl comment form, which is how these pages render them;
`.github/workflows/Scripts.yml` runs every one to completion, and every
published build of this site runs them again:

```bash
julia --project=verification -e 'using Pkg; Pkg.instantiate()'   # once, Julia ≥ 1.11
julia --project=verification verification/landau-damping-1d1v.jl
julia --project=verification verification/plasma-oscillations-1d1v.jl
julia --project=verification verification/two-stream.jl
julia --project=verification verification/bump-on-tail.jl
julia --project=verification verification/plasma-echo.jl
julia --project=verification verification/bgk-equilibrium.jl
julia --project=verification verification/collisional-damping.jl
julia --project=verification verification/wakefield.jl
julia --project=verification verification/scheme-comparison.jl
```

* `landau-damping-1d1v.jl`: linear Landau damping, the recurrence of its modes,
  and the strong case at α = 0.5, on uniform and non-uniform velocity grids.
* `plasma-oscillations-1d1v.jl`: the Bohm–Gross frequency, and the energy to
  t = 3000 on both grids.
* `two-stream.jl`: the growth rate across the unstable branch, against the warm
  kinetic root.
* `bump-on-tail.jl`: Arber and Vann's beam, through saturation by trapping.
* `plasma-echo.jl`: the echo with the field off and on, and what four schemes
  keep of it.
* `bgk-equilibrium.jl`: a nonlinear stationary state, and where the solver
  moves it.
* `collisional-damping.jl`: `BGK` in the loop, against its own dispersion
  relation.
* `wakefield.jl`: a laser driving a plasma wave, against linear theory.
* `scheme-comparison.jl`: the advection schemes ranked on the Landau rate and on
  positivity. It is advisory and exits 0.

Trapping, reversibility, the Doppler shift, convergence and the invariants are
runs of the test suite alone, with no study of their own.

## What the tests assert

The studies' headline claims are asserted by the test suite rather than left in prose.
CI runs this on every pull request; locally it is behind an environment variable
so that a default `Pkg.test()` stays instant:

```bash
VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

The analytic side of those claims is computed rather than quoted:
[`test/dispersion.jl`](https://github.com/korzhimanov/Vasilek.jl/blob/master/test/dispersion.jl) solves the kinetic dispersion relation through the plasma
dispersion function, which gives the Landau roots at any `k` — they used to be
three constants typed into the test file — the growth rate of *warm*
counter-streaming beams, which is what the runs contain, and the growing root of
the bump-on-tail distribution, on the continuum and on the grid's own
centred-difference field. The cold two-stream closed form remains as its
zero-temperature limit and is checked as such. With the `BGK` operator in the
loop, the same `Z` gives the collisional relation: a 3×3 closure on the three
moments the operator restores, checked against the Landau root at `ν = 0`, the
Krook model's closed form, and the Chapman–Enskog fluid at large `ν`.
`Z(ζ) = i√π·erfcx(−iζ)` is written through the entire Faddeeva function, so one
expression covers the whole plane, and `test/test_dispersion.jl` holds it to
identities rather than to stored numbers: `Z` against its own Cauchy integral
where that converges, to 3.5e-15, and the susceptibility against a trapezoid
that never calls `Z`, to 1.3e-14.

Every row below is asserted in `test/test_verification.jl`, except the field-off
echo, in `test/VlasovSolver/test_echo.jl`, and the collisional relation, in
`test/test_dispersion.jl`.

| claim | asserted | measured |
|---|---|---|
| Landau damping rate, k = 0.3, 0.4, 0.5 | γ within 3% of the root of the kinetic dispersion relation | 0.70%, 1.14%, 0.95% |
| Landau real frequency, same three k | ω within 1% of the same root | 0.08%, 0.22%, 0.26% |
| and the runs are linear enough for that to be the right target | bounce phase under 2 in every fitting window | 1.70, 0.45, 0.21 |
| the damping fit does not depend on its window | two windows agree within 3% | 0.17%, 0.05%, 0.04% |
| nor does the frequency fit | two windows agree within 1% | 0.066%, 0.000%, 0.300% |
| the Vlasov–Poisson flow runs backwards | round-trip error after v → −v falls by at least 3.5× when the grid is halved, at α = 0.05 | 2.82e-3 → 6.48e-4, ×4.3: second order, the limiter clipping the peak |
| unless the flow has out-run the grid | at α = 0.5 the same refinement recovers nothing | 0.167 → 0.157, ×1.06 |
| and reversibility is bought with positivity | the two most reversible schemes are the two that drive f negative | 0.090 and 0.108 at f = −0.093 and −0.057, against PFC's 0.167 at f ≥ 0 |
| strong Landau damping, α = 0.5 | γ₁ within 0.27…0.30, around the literature's 0.281…0.292, which it quotes as −γ; γ₂ within 0.070…0.090, just under the cited 0.0815…0.0858, and refining moves it toward them | 0.2863, 0.0787 (0.0813 at twice the resolution) |
| and the non-uniform velocity grid agrees | both rates within 2% of the uniform grid at the same Δt | 0.06%, 0.5% |
| a drifting plasma damps the same way | boosted mode matches the rest-frame one once the Doppler phase is removed at each sample's own time, and the fitted rates agree within 0.1% | 1.7e-3, 1.8e-3 rad, 0.014% |
| collisional Landau damping, `BGK` at ν = 0.1, 0.3, 1 | γ within 1.5% and ω within 0.3% of the collisional root on the grid's field, and γ *falling* as ν rises | +0.37%, +0.49%, +0.76%; +0.05%, +0.08%, +0.03%; 0.1547 → 0.1341 → 0.1074 → 0.0660 at ν = 0, 0.1, 0.3, 1 |
| which is BGK's root and no other operator's | a relaxation restoring only n, or only n and u, lands on its own root and over 1.5× BGK's rate | 0.31%, 0.23% from their own; 3.3×, 2.1× BGK's |
| and the mode energy conservation adds | a matrix pencil of the field finds a real exponent within 0.2% of `heat_mode_root` at ν = 1 | 0.43377 against 0.43369 |
| the collisional relation is the right one | the Landau root at ν = 0, Krook's closed form for density alone, and Chapman–Enskog as ν → ∞ with gaps closing as 1/ν² | 2.3e-14, 4.7e-16; gaps ×3.85–4.08 per doubling |
| and `BGK` holds the energy on any window | at ν = 1 the energy holds to 1e-6 on ±8 and on ±4, where the sampled Maxwellian of 0.1 cools | 8.0e-9, 8.0e-8; −5.3% by t = 60 |
| trapping stops the damping on the bounce time | ω_B·t₀ between 6.5 and 8.5 at four amplitudes, and t₀ ∝ α^(−1/2) | 7.09–8.00, slope −0.556 |
| each mode recurs at its own 2π/(kΔv) | within a plasma period, for the seeded mode and the harmonic it generates | 128.6 vs 125.7, 64.3 vs 62.8 |
| a plasma echo, field off, is the closed form's | pointwise within 0.3% of the peak, the peak within two steps of its own, the sign reversing with the kick | 0.12%, t = 15.28 on both |
| out of a mode no moment could see | the echo over 10⁴ times the seeded mode's density at the kick | 2.4·10⁴ |
| and what a scheme keeps of a filament is what it returns | at 64 × 121, SemiLagrangian < PFC < LaxWendroff < Upwind; upwind returns under 70% | 0.11%, 1.03%, 4.5%, 45%; 55% |
| the loss is truncation | PFC's error falls at least 5× per halving of Δv | 7.4×, 6.9× |
| a plasma echo, field on, is second-order kinetic theory's | pointwise within 1% of `echo_second_order`, the peak within 1% of its size and two steps of its time | 0.48%, 0.48%, t = 29.05 on both |
| which is not the field-off echo | the closed form's peak over 1/0.6 times the run's and 0.5 later | 2.1×, 1.14 |
| and the residual is the run's | it falls 2.5× from Nx = 64 to 128 | 3.3× |
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
| two-stream at the cold boundary, kv₀ = 1.0 | γ within 8% of the warm root, where the cold form gives exactly zero | 3.16% |
| two-stream stability boundary, kv₀ = 1.2, 1.6 | no growth, where cold and warm theory both give γ = 0 | decays to 0.052, 0.000 |
| bump-on-tail growth, Arber and Vann's beam | γ and ω within 0.05% and 0.01% of the kinetic root on the grid's field, 0.5% and 0.05% of the continuum's; the wave travels with the beam | 0.016%, 0.003%; 0.131%, 0.014% |
| it saturates by trapping, whatever the seed | ω_B/γ between 1.8 and 2.1 at the peak; a seed 1000× larger peaks within 0.5%, ln(1000)/γ earlier | 1.951; 0.04%, 35.00 vs 34.91 |
| and the trapped beam carries the field | the amplitude swings with a period of 1 to 1.6 bounce periods; the flank's slope in ⟨f⟩ falls at least 3× | 1.316; 20× |
| plasma oscillations, uniform grid | abs(Δε/ε) < 0.5% at t = 3000 | 0.38% |
| plasma oscillations, non-uniform grid | abs(Δε/ε) < 6% at t = 3000 | 4.75% |
| and the energy is one instant's | on the uniform grid ε moves by under 3e-4 of itself within a plasma period | 6.6e-5, where summing the kinetic energy after the kick swung it by 1.2e-3 |
| the laser wakefield is a plasma wave | ω within 2% of Bohm–Gross √(ωₚ² + 3Tk²) | 0.23% |
| the driver travels at the grid's group velocity | pulse speed within 2% of `vg_pulse`, at two resolutions | 0.73%, 1.51% |
| and refining moves it toward the continuum | abs(v − √(1−n)) falls when Δx is halved | 0.9487 − 0.8858 → 0.9487 − 0.9261 |
| its wavelength is the driver's | λ within 3% of 2π√(v² − 3T)/ωₚ, for the measured *and* the predicted driver | 0.37%, 1.12% |
| it is phase-locked to the pulse | ω/k within 2% of the pulse velocity, measured independently | 0.58% |
| its size is the one linear theory gives | peak within 10% of `linear_wake`, pointwise rms under 15% | 6.4%, 6.4% |
| and it is the laser that made it | amplitude ∝ a₀² within 8%; ≥10× the unlit control; ≥4× behind the pulse over ahead | 3.3%, 36×, 9.4× |

## Measurements by study

Each section says what is claimed, what measures it, and the numbers. The
measurement journals these were gathered from, as they stood through commit
`f5032f9`, are in the git history of `CHANGELOG.md`.

### Landau damping and dispersion

One `k` with one fitted `γ` is a single point on a curve, and cannot separate a
solver that reproduces the dispersion relation from one that lands near a value
at one wavenumber. Both roots are measured at k = 0.3, 0.4 and 0.5 against
`landau_root`, which solves the kinetic dispersion relation rather than quoting
it (testset "Landau damping and dispersion"; `damping_rate`,
`oscillation_frequency` and `trapping_phase` in the harness). The real frequency
is four to nine times the sharper of the two, since it comes from counting nulls
where `γ` comes from fitting an amplitude that numerical dissipation also acts
on, so it is held to 1% where `γ` is held to 3%.

Three constraints shape the runs, and each cost time to rediscover:

* **The velocity window must contain the resonance** at `v = ω/k` (2.83, 3.21,
  3.87), since that is where the damping comes from; the window then fixes Δt,
  because the fastest row runs at `max|v|·Δt/Δx` against the scheme's Courant
  limit. The k = 0.5 case runs at 0.815.
* **The fitting window stops before the floor** where recurrence and round-off
  take over, and the floor arrives earlier in e-foldings the weaker the damping,
  which is why each k has its own window; at α = 1e-2, where k = 0.3 had decayed
  only threefold by t = 50, fitting past it read γ 22% low.
* **The runs are linear.** At α = 1e-2 the k = 0.3 case is not: trapping stops
  its damping inside the fitting window (next section), and the γ it reported sat
  *below* the root, two errors cancelling, where dissipation can only push it
  high. At α = 1e-3 the rate is flat to t = 95; the fitting window is [10, 90],
  thirty maxima, and the bounce phase `trapping_phase` reports over each window,
  1.70, 0.45 and 0.21 at k = 0.3, 0.4 and 0.5, is asserted under 2.

`damping_rate` fits a line through the maxima of `log ε_e` only. Fitted through
every sample, the `log cos²` nulls pull the slope by an amount that depends on
where the window starts; that is what once read 0.1498 at k = 0.5. Both
estimators are held to their windows: two windows agree in `γ` to 0.17%, 0.05%
and 0.04%, and in `ω` to 0.066%, 0.000% and 0.300%. `ω`'s 1% sits in a gap.
Below them is the estimator's own floor: a null is located only to within Δt,
so two windows disagree by about `Δt/span` whatever the physics does, 0.4% at
k = 0.5. Above them is the failure tested for: losing one null of ten rescales
the spacing by 11%.

The study renders the k = 0.5 run's spectrum, which peaks at 1.4137 against the
root's 1.4157.

### Trapping and the bounce time

The linear cases are linear because their amplitude was chosen to make them so.
At larger amplitudes trapped particles stop the damping, which is O'Neil's
result, and reproducing *when* it happens is a statement about the nonlinear
term that no damped-mode test can make (testset "Trapping stops the damping, on
the bounce time", at k = 0.3). At α = 1e-2 the local rate between maxima two
apart holds near 0.013 to t ≈ 30, reaches zero at t ≈ 73 and goes negative
after it, the field growing again as the trapped population sloshes; at
α = 1e-3 the same rate stays between 0.0121 and 0.0129 to t = 95.

With `ω_B = √(kE₀) = √α` here, the arrest time should go as α^(−1/2), and the
bounce phase `ω_B·t₀` at which it happens should not move. The arrest time is
interpolated between the last positive local rate and the first non-positive
one, so it is not quantised to the spacing of the maxima:

| α | t₀ | √α·t₀ |
|---|---|---|
| 0.005 | 113.2 | 8.00 |
| 0.01 | 73.4 | 7.34 |
| 0.02 | 50.5 | 7.14 |
| 0.04 | 35.5 | 7.09 |

A factor of eight in amplitude moves the arrest by 3.2, where α^(−1/2) predicts
2.83; the fitted slope of log t₀ against log α is −0.556. The same scan
calibrates `trapping_phase`, the guard the linear cases are held under: at the
arrest times it reads 4.26 to 5.72, the damping departs from the root at about
2.8, and the linear runs sit at 0.21 to 1.70.

### Recurrence

On a discrete velocity grid a phase-mixed mode does not stay mixed: the sheared
modes rephase, and the mode returns at `T_R = 2π/(kΔv)`, a time set by Δv alone
that no scheme can move. In the Landau run the seeded k = 0.5 mode returns at
128.6 against 125.7, at 52% of its initial amplitude, and the k = 1 harmonic it
generates returns at 64.3 against 62.8, from a round-off floor of 6.25e-17; over
t ∈ [58, 72] the harmonic is 21 times the seeded mode (testset "Each mode recurs
at its own 2π/(kΔv)"; the driver records any modes asked of it through
`modes`). What the Landau study's figure shows near t = 63 is therefore the
harmonic's recurrence, not the fundamental's.

With the field off the Vlasov equation has an exact solution, and the same
statement is made in `test/VlasovSolver/test_free_streaming.jl`: a modulated
Maxwellian's density mode decays as `exp(−k²t²/2)`, a Gaussian rather than an
exponential, because the mode is sheared into velocity structure rather than
dissipated. Over t ≤ 4 the mode amplitude is off by 2.7e-6 for cubic
`SemiLagrangian`, 3.3e-5 for PFC, 1.2e-3 for `LaxWendroff` and 1.7e-2 for
`Upwind`, and the ordering is asserted. On a coarse velocity grid, Δv = 0.4 and
`T_R = 31.4159`, cubic `SemiLagrangian` and PFC peak at t = 31.42, recovering
1.0000 and 0.9991 of the initial amplitude. The Courant limit binds here as it does not in
one-dimensional advection: the fastest velocity row carries `max|v|·Δt/Δx`, and
everything runs at 0.611.

### The plasma echo

The echo is the one place a kinetic theory beyond linear order is held to a
number. A seed at k₁ = 1, a velocity kick at k₂ = 3/2 when t = 5, and around
t = 15 a mode nobody seeded appears at k₃ = 1/2: phase mixing is reversible,
and the echo is made of information no moment of `f` can see. So it measures
what a scheme keeps of filaments finer than any it can show, which free
streaming, measuring only the decay, never did.

**Field off** (`test/VlasovSolver/test_echo.jl`, fast suite). The echo has a
closed form, exact in both amplitudes,
`−iα·J₁(k₃ε(t − τ))·exp(−(k₃t − k₂τ)²/2)` (`echo_closed_form` in
`test/echo.jl`). The Bessel function is not a refinement: its small-argument
limit is 14.6% high at these amplitudes. The sign is asserted separately, since
the direction of the kick changes it and nothing about the peak. PFC at
128 × 241 follows the closed form to 1.25e-3 of the peak, pointwise, out of a
seeded mode the kick found at 1.84e-6, 2.4e4 times smaller than the echo it
turns into. The loss converges at orders 2.89 and 2.78 as Δv halves, and, at
64 × 121, it ranks the schemes: cubic `SemiLagrangian` 1.07e-3, PFC 1.03e-2,
`LaxWendroff` 4.48e-2, and `Upwind` 0.449, which returns 55% of the echo and
loses most of the rest in the x-sweep rather than in the kick.

**Field on** (testset "The plasma echo, with the plasma's own field", extended
suite). `echo_second_order` composes three linear responses of the Maxwellian:
the seed's filament screened at k₁, the kick screened at k₂ together with the
lifetime of the field it induces, and the echo's own density polarising the
plasma at k₃, solved as a Volterra equation. Its pieces are checked before any
run is held to it: the resonant dielectric function is `dielectric` to 1.2e-16,
the filament's closed form agrees with a time-domain solve to 4.9e-6, and with
the field off the theory is the closed form's small-ε limit to 1.1e-14. The
self-consistent run (`self_consistent_echo`, 128 × 241, α = ε = 0.01, τ = 10)
follows it to 4.8e-3 of the peak, pointwise, the ringing after the peak
included, with the peak 0.48% low on the same sample, t = 29.05. The field-off
form would put that peak 2.1 times higher and 1.14 later. The residual is the
x-sweep's: 1.6e-2 at Nx = 64, where halving Δv or Δt instead leaves it at
1.5e-2 and 1.7e-2.

### The BGK equilibrium

The equilibrium study is the converse of every other run: it starts on a
nonlinear stationary state and asks how little the solver moves it. The state
is a function of the particle energy in the potential `−ψ cos kx`, on an ion
background built to hold it (`bgk_equilibrium` in the harness; the driver keeps
the ions it is handed as `nᵢ`). At ψ = 0.5 two thirds of the particles are
trapped, so this is also the first run whose physics they carry (testset "A
nonlinear equilibrium stays put, and gives on the separatrix").

The Maxwell–Boltzmann equilibrium, analytic across the separatrix, holds `f` to
1.52e-3 of its peak and the field to 3.85e-3 through t = 50. That is a steady
drift, not an oscillation about the equilibrium, and it is the scheme's
dissipation: L² falls by 1.13e-3 and the entropy rises by 1.27e-4, where the
exact flow keeps both. The error spreads over phase space and converges at
PFC's third order, 7.6 times less in `f` and 6.1 in the field per halving;
halving Δt alone moves neither.

Give the trapped particles a temperature of their own and `F` keeps its value
but not its slope across the separatrix. The worst cell is then *on* the
separatrix at both resolutions, where the smooth case's is 13 and 24 cells away,
the error is twelve times the smooth one, and it converges at first order, 1.75
per halving. Trapped-particle structure lives or dies there.

The teeth: ions built for a potential 10% off put the field 49% and 59% off at
once, and `f` drifts 24 and 27 times further; ions built on the Poisson sign the
normalization notes once gave hold the reversed field, `E/E₀ = −1.000`, and the
run leaves the equilibrium by 16% of its peak.

A trapped population colder than the bulk peaks at 1.79, above the f = 1 the
harness once bounded PFC by. Bounded by `f₀`, as every run now is, it holds to
0.99% of its peak. At the old bound the limiter, handed data above its `fmax`,
turned inside out and took `f` 43.6% of its peak away, without an error.
`PFCNonUniform` now checks its data as PFC does and refuses that bound on the
first call; built unchecked, it reproduces the 43.6%, and the extended suite
asserts both. Over the two-stream runs carried into saturation and the three
equilibria, 3.8 million calls of the checked schemes stay inside their bounds,
sub-steps included.

### Reversibility

Vlasov–Poisson is invariant under `(t, v) → (−t, −v)` with `E` unchanged. Run
forward, flip the velocity axis, run forward again, flip back: an exact scheme
returns the initial state, and what a real one loses is its own dissipation
(testset "The Vlasov-Poisson flow is reversible, and what breaks it";
`vlasov_poisson` returns `f` so that the question can be asked of the
distribution rather than of a moment of it).

* At α = 0.05 the round-trip error is 2.82e-3, 6.48e-4 and 1.58e-4 as the grid
  halves, ratios 4.3 and 4.1: second order, so the irreversibility is the
  scheme's and converges away. Bounded at the initial maximum, PFC's limiter
  clips the peak cell, which costs the third order it shows when bounded at 1
  (ratios 5.7 and 7.1); the test asks for 3.5.
* At α = 0.5 it is 0.167 and stays 0.157 when the grid is halved, a factor of
  1.06. By then the flow has folded the distribution into filaments finer than
  Δv, and the information needed to run the film backwards is not on the mesh to
  be refined. Every nonlinear run in the suite is in that regime by the time it
  is interesting.
* Ranked by round trip at α = 0.5: `LaxWendroff` 0.090, cubic `SemiLagrangian`
  0.108, PFC 0.167, `Upwind` 0.332. The first two get there by driving `f` to
  −0.093 and −0.057 against a peak of 0.6, over the round trip's two runs of
  t = 20 on 64 × 121 cells, where PFC stays at +5e-10 and
  `Upwind` pays a fifth of its L² norm instead. Reversibility and positivity are
  the same trade-off from two sides, reached through a quantity that has
  nothing to do with the damping rate the scheme comparison ranks them by.

The same statement with the field off, a rigid rotation run forward, flipped,
run forward and flipped back (`test/VlasovSolver/test_strang_splitting.jl`),
separates the schemes far more sharply than the forward error does. At
N = 64 → 128: cubic `SemiLagrangian` 1.01e-3 → 1.23e-4, `LaxWendroff`
3.77e-3 → 4.94e-4, PFC 1.95e-2 → 3.61e-3, and `Upwind` 4.04e-1 → 2.56e-1.
`Upwind` is there because it fails: 40% of the peak, smeared past recovery, and
refining barely helps because its dissipation is first order.

### Strong Landau damping

The α = 0.5 case is the standard nonlinear benchmark, and it is quoted by two
numbers: Cheng and Knorr (1976) give γ₁ = −0.281 and γ₂ = 0.084, later work
−0.292 with 0.0815 and −0.2918 with 0.08584, the literature quoting the decay
with a negative sign. Here γ₁ is the rate of decay, held within 0.27…0.30, and
γ₂ the rate of regrowth, within 0.070…0.090 (testset "Strong Landau damping:
the damping stops and reverses"). Measured on 128 × 241 with Δt = 0.025:
γ₁ = 0.2863 over the four maxima of the decay proper, γ₂ = 0.0787 over the eight
of the regrowth. Both windows are conventions, and the testset says so with
numbers: γ₁ reads 0.3793 over three maxima and 0.2283 over five, because the
envelope steepens and then flattens into the trapping plateau, so a straight
line through it depends on how much of the curve is inside the window.

Refinement is what makes the agreement more than a coincidence of one grid: γ₂
goes 0.0714 → 0.0787 → 0.0813 as Δx, Δv and Δt halve, toward the published value
rather than away from it. The non-uniform velocity grid, the only path the suite
has toward an adaptive mesh, reproduces both rates to 0.06% and 0.5%, on a run
that reaches the sharp gradients its limiter exists for.

### The Doppler shift

Every other Vlasov–Poisson case starts from a distribution symmetric in `v`, so
the mean velocity is zero throughout and the drifting half of the solver is
never exercised. A boost by `u` carries `f(x, v, t) → f(x − ut, v − u, t)` and
`E(x, t) → E(x − ut, t)`, so the mode's history should be the same complex
function times `exp(−ikut)` (testset "A drifting plasma damps the same way,
Doppler-shifted", at u = 0.5). Nothing in the discretisation is Galilean
invariant, since the grid does not move and the boosted Maxwellian sits on it
asymmetrically, so the agreement is a measurement: over the eleven maxima in the
fitting window, amplitudes within 1.7e-3 and phases within 1.8e-3 rad, with
fitted rates 0.014% apart and frequencies identical to the estimator's
resolution. Without the Doppler correction the phase differs by up to 2.88 rad,
which is what says the correction is doing work.

The phase is compared at each sample's own time, `t[k] + Δt/2`, where the
driver records the field ([driver notes](driver-notes.md)); removing the Doppler
factor at `t[k]` instead reads 8.0e-3, of which 6.25e-3 is that half-step and
not the grid. The rates agree to 0.014% because the energies are the cell-width
sums the schemes conserve; the trapezoid, which books a travelling wave's passage
across the seam as energy, put them 0.136% apart. The comparison is made at the
maxima rather than pointwise: both `|E_k|` and `ε_e` pass through deep nulls,
and two nulls landing a time step apart read as a 900% discrepancy.

### Convergence under refinement

Agreement at one resolution can be a coincidence of that grid; the testset
"Landau damping converges under refinement" asks that it improve.

In the Landau run, halving Δx, Δv and Δt together, so that the Courant number
stays at 0.81 and only the discretisation moves, gives γ errors of 6.47%, 1.33%,
0.71% and 0.47% for Nx = 32 to 256. The L² dissipation over the same ladder
falls by 11.0, 7.9 and 7.8 times against the 8 a third-order scheme predicts,
which identifies the residual error in γ: the fitted rate is the physical
damping plus the scheme's own, which is why every measurement sits on the high
side of the root rather than scattering about it. The testset runs the first
three rungs and asks that the error fall at each and end under 1%, and that the
dissipation fall between 4 and 16 times per halving.

The splitting is measured with a real scheme against an analytic answer in
`test/VlasovSolver/test_strang_splitting.jl`, on a rigid rotation, refining Δx,
Δv and Δt together over N = 32, 64, 128: cubic `SemiLagrangian` goes
6.37e-3 → 5.72e-5 at orders 3.60 and 3.19, PFC converges at 2.26 and 2.22, and
`LaxWendroff` at 1.84 and 1.96. Each scheme keeps its own spatial order rather
than being dragged to the splitting's second. That is not the problem being
easy: a rigid rotation factors into shears, which is what the splitting does, so
the commutator error is unusually small here and the scheme is what is left.

### Invariants

Over the k = 0.5 Landau case, 875 steps, mass drifts 1.6e-15 and momentum stays
at 2.2e-15 on a mass of 25.1: round-off, both exact conservation laws the
discrete scheme also satisfies. L² falls 1.0e-5 and the entropy rises 7.5e-6,
monotonically at every step. Those two are *not* conserved and are not asserted
as if they were: an exact Vlasov flow preserves both, and the drift is numerical
dissipation, so what is asserted is the direction (testset "The invariants of
the split system"; `test/test_invariants.jl` for the kernels alone).

The invariants and the energies are the cell-width sums `Σ f ΔvΔx` and
`Σ e² Δx`, the quadrature a flux form conserves, not the trapezoid
([driver notes](driver-notes.md)). On the same run the trapezoid reports 2.4e-4
of mass drift and 7.3e-4 of momentum, eleven orders of magnitude of apparent
non-conservation that belongs to its endpoint weights. On the bump-on-tail run,
whose wave travels across the seam, the trapezoid's `ε`, 61.8, swung by −0.346
to +0.175 over t ∈ [60, 100] on 64 × 361 cells and by −0.185 to +0.093 at
Nx = 128, first order in Δx, where the sums drift by 9.7e-3 and 1.0e-3.

### Collisional damping

The collisional study is the first run with the collision operator in it. `BGK`
had been checked on a single velocity line, and never with a field. Put between
the kicks at three rates, it lands on its own dispersion relation, which says
something one might not guess: collisions that restore density, momentum and
energy *weaken* the damping, from Landau's 0.153 towards a fluid wave that damps
only by conducting heat. A relaxation that restored less would damp harder, and
the testset ("Collisional Landau damping: BGK in the loop") runs two such
relaxations to show that each lands on its own root.

**The relation.** Linearised, BGK relaxes a mode towards
`F₀[n₁ + u₁v + T₁(v² − 1)/2]`, and the three moments close the kinetic equation
into a 3×3 determinant built from `Z` (`collisional_determinant`,
`collisional_root` and `heat_mode_root` in `test/dispersion.jl`). At `ν = 0` the
determinant *is* the collisionless dielectric function, to 2.3e-14; restricted
to the density it is the Krook model's closed form, to 4.7e-16. At large `ν` it
approaches the Chapman–Enskog fluid, derived separately: no viscosity in one
dimension, heat flux `−(3nT/ν)∂ₓT`, and a cubic whose roots are the Langmuir
pair and a heat-conduction mode. The gaps in ω, γ and the heat rate close as
1/ν², 3.85 to 4.08 times per doubling, as the next order of the expansion says
they should. Where the root sits above the real axis, at `ν ≥ 0.3`, a trapezoid
over the mode's own moments, which never calls `Z`, confirms it (2.3e-14). At
k = 0.5 the root's γ goes 0.1534 → 0.1327 → 0.1062 → 0.0652 at ν = 0, 0.1, 0.3,
1, towards a wave damped by conduction alone at `3k⁴/(ν(1 + 3k²))` and
oscillating at the adiabatic `√(1 + 3k²)`. Restoring fewer moments strengthens
the damping: dγ/dν at ν = 0 is −0.243 for BGK, +0.114 without energy and
+0.684 for Krook. The collisional roots take a looser secant stop than the
collisionless ones, 1e-9 rather than 1e-12: `J₃` and `J₄` come from `1 + ζJ₂`,
which cancels when ν is large, and against quadrature they carry 1.2e-10 of
their size at ν = 20. `heat_mode_root` searches the imaginary axis only above
`Im ζ = −2`; further down the continuation of `Z` grows as `exp(|ζ|²/2)`.

**BGK.** On 64 × 161 over v ∈ [−8, 8] with Δt = 0.04, each run lands on the
root on the grid's field (`collisional_landau` in the harness): γ reads +0.37%,
+0.49% and +0.76% above it at ν = 0.1, 0.3 and 1, ω +0.05%, +0.08% and +0.03%,
and the damping falls as ν rises. The residue in γ is the x-sweep's
dissipation: halving Δx and Δt at ν = 1 takes it to 0.03%, while halving Δv
moves γ by 0.01%.

**Partial relaxations.** `PartialBGK` relaxes at the background temperature,
restoring only the density (the Krook model) or the density and the drift. Run
in place of BGK, each lands on the root `collisional_root` gives for what it
restores, to 0.31% and 0.23% in γ, at 3.3 and 2.1 times BGK's rate. So the
relation is checked for three operators, and a BGK that lost its temperature or
its drift would fail by a factor.

**The heat mode.** Energy conservation gives the plasma a third mode at each
`k`, a temperature perturbation that does not oscillate. The density seed
excites it at 21% of the wave's amplitude at t = 4, and it decays six times
faster than the wave, so the maxima and nulls never see it. A matrix pencil of
the field's `k` mode, `mode_exponents`, returns the wave pair and a real
exponent, 0.43377 against `heat_mode_root`'s 0.43369 at ν = 1; its wave agrees
with the maxima and nulls to 0.06% in γ and 0.04% in ω. The pencil takes the
*left* singular vectors: the right ones return every exponent conjugated, which
on a test signal swapped a travelling wave's direction and nothing else.

**The window.** The operator relaxes each line to its *discrete* Maxwellian, the
`exp(a + bv + cv²)` whose density, momentum and energy on the grid equal the
line's, so it conserves all three to round-off on any window. At ν = 1 the
energy holds to 8.0e-9 on ±8 and to 8.0e-8 on ±4, and the entropy rises
monotonically. The sampled Maxwellian it used before, still there as
`BGK(τ; conservative = false)`, took the temperature by the trapezoid over the
window and put back a Maxwellian over the whole line, so a window that cut the
tail cooled the plasma at every relaxation, by 1.2e-3 of the energy: by t = 60
on ±4 the run has lost 5.3% of its energy, T has fallen to 0.949, and γ reads
3.0% low.

### Plasma oscillation frequency and energy

A solver oscillating at entirely the wrong rate passes an energy test without
difficulty, so the frequency is measured as well (testsets "Plasma
oscillations, Bohm–Gross frequency" and "Plasma oscillations, energy
conservation"). At k = 2π/100, ω = 1.005719 against `√(1 + 3k²)` = 1.005904,
0.018%, where the cold `ωₚ = 1` is 0.57% away. The factor of thirty is the
point: the test fails if the thermal correction is dropped rather than merely
preferring it there. That costs a second at t ≤ 200, where the energy test needs
t = 3000 and forty times as long. Two windows agree to 0.006%.

The study renders the spectrum of the mode the driver records, projected on the
standing wave's phase and zero-padded; it peaks at 1.00528 against Bohm–Gross's
1.00590.

Over t ≤ 3000 the energy drifts 0.38% on the uniform velocity grid and 4.75% on
the non-uniform one, against 0.5% and 6%. The non-uniform figure was about 12%
while `PFCNonUniform` took its limiter from the grid's global spacing ratio
rather than each cell's own. Within a plasma period `ε` moves by 6.6e-5 of
itself, asserted under 3e-4, because its kinetic energy is the mean of the two
sides of the kick and so is sampled at `t[k] + Δt/2` with the field
([driver notes](driver-notes.md)). Summed after the kick it swung by 1.2e-3,
half the kick's work counted early, and the agreement at t = 3000 was then a
point on that swing landing near the drift.

### The two-stream instability

Every other case here is a damped or neutral mode. A growth rate catches a class
of error damping cannot: a sign flip in the field push turns damping into growth
and growth into damping, so a suite of damped cases alone is half-blind to it
(testset "Two-stream instability"; `two_stream` and `growth_rate` in the
harness, `two_stream_warm` in `test/dispersion.jl`).

The beams are Maxwellians at ±v₀ with `vt = 0.3`, and the reference is the warm
kinetic root. The cold relation, `1 = ½/(ω − kv₀)² + ½/(ω + kv₀)²`, is a quadratic
in ω², `u± = [(2a² + 1) ± √(8a² + 1)]/2` with `a = kv₀`, unstable exactly when
`a < 1`, with `γ = √(−u₋)`; it stays as the warm root's zero-temperature limit,
checked against its relation to 1.4e-14. At a = 0.4, 0.6 and 0.8 the runs read
−0.63%, −1.95% and +0.09% from the warm root; from the cold one they read
−2.10%, −3.15% and +0.30%, half a tolerance spent on a known approximation.
`γ(a)` is not monotone and peaks near `a = √(3/8)`, so reproducing all three is
a statement about the branch: a solver that merely amplified what it was given
could not put the maximum in the right place.

Finite temperature can make a beam grow faster than a cold one: the warm root
crosses above the cold between a = 0.75 and 0.8, and past a = 1 the cold branch
is identically zero while warm beams still grow. At a = 1.0 the warm root
predicts 0.09823 and the run gives 0.09513, 3.16%, held to 8% because the ripple
below is worst where γ is smallest: over the band's upper edge `hi` ∈ {1, 2, 3,
5} the fit moves between −6.05% and −3.16%. At a = 1.2 and 1.6, where both
theories give γ = 0, the mode decays to 0.052 and 0.000 of its initial energy.
Doubling the beams' temperature at a = 0.6 lowers γ by 4.4%, where warm theory
predicts 4.4%, and each run is within 3% of its own root. The test asserts that
relative statement, two temperatures at one wavenumber, rather than any sign at
a single wavenumber.

`growth_rate` fits every sample in a band of amplitude rather than of time: the
unstable root is purely imaginary, so the mode grows without oscillating and
there are no maxima to fit. That makes it transferable across the branch: at
a = 0.4 fixed time windows give 10.65%, 6.43%, 3.99% and 0.19% from the cold root
depending on where they are put. `ε_e` is not one exponential either. The
cold quadratic's oscillating pair, `±ω₊`, puts a ripple on it whose size relative
to the growing mode falls only as `exp(−γt)`: at a = 0.6 the instantaneous rate
oscillates with period 4.5 against `2π/ω₊` = 4.626, between 0.21 and 0.41 about a
cold γ of 0.353. Fitting over two beat periods rather than one cuts the spread
over start points from 41.4%, 14.4% and 22.3% to 9.4%, 5.6% and 4.8%, and the
amplitude band already spans 1.15 to 2.22 periods. Under the band, cooling the
beams at a = 0.6 from `vt = 0.6` to 0.15 brings the rate to the cold one
monotonically from below, −7.44% to −2.10%. An exploratory fixed window,
t ∈ [10, 20], had it cross the cold value, +0.40% and +0.56% at `vt` = 0.2 and
0.15: changing `vt` changes γ slightly, which moves the ripple's phase within a
fixed window and drags the fitted rate with it.

The box holds 96 cells by construction, `range(Δx; step = Δx, length = Nx)`: a
float range built from its endpoints came up a point short at a = 0.4 and gave
the box another fundamental, a = 0.4042. Past their fits the runs grow until the
velocity sweep would cross its Courant limit, where the driver splits the step,
so they saturate rather than diverge, `ε_e` peaking at 106 for a = 0.6. Every fit
ends at a velocity Courant number of 0.53 to 0.73, before the first split.

### The bump-on-tail instability

The bump-on-tail study is the first run held to numbers past its linear phase.
The beam is Arber and Vann's [J. Comput. Phys. 180, 339 (2002)],
`F = [0.9·exp(−v²/2) + 0.2·exp(−2(v − 4.5)²)]/√(2π)`, over one wavelength of
k = 0.3: a Langmuir wave travelling with a weak beam over a warm bulk, where the
two-stream case grows as a standing wave. At this beam's temperature the
instability is mostly reactive; the temperature cuts the rate by 15% from the
cold beam's, and the run has to resolve that cut. Measured on 64 × 361 with
Δt = 0.05 (`bump_on_tail` and `mode_rates` in the harness; testset
"Bump-on-tail: the beam grows a wave, and the wave traps the beam"):

* **Growth.** γ = 0.19784 and ω = 1.00108 over t ∈ [30, 50], from a seed of
  1e-6: 0.131% and 0.014% below the kinetic root, 1.001218 + 0.198098i. Nearly
  all of the gap is the Poisson solve's centred difference, which returns
  `sin(kΔx)/(kΔx)` of a mode's field. Put that factor on the susceptibilities,
  as `dielectric` does given `Δx`, and the root sits 0.016% and 0.003% from the
  run. The rate is set in x, not v: Δv from 0.025 to 0.2 moves γ by 0.017% at
  most, while Nx = 32 puts it 0.69% off the continuum and Nx = 128 0.020%. The
  phase of the mode amplitude shows the wave travelling with the beam.
* **Saturation.** The mode grows until it traps the beam feeding it, and peaks
  at `|E_k|` = 0.4979 at t = 74.32, where the bounce frequency is 1.951γ
  (1.945–1.952 over the grids tried). A seed a thousand times larger peaks at
  0.4981, 35.00 earlier, against ln(1000)/γ = 34.91.
* **Trapping.** The trapped beam swings round the well and the field swings with
  it: past the peak the amplitude falls to 0.2645 and comes back to 0.4137, with
  a period of 21.40, 1.32 bounce periods at the bottom of the well. The averaged
  distribution loses the slope the growth ran on: the steepest slope of the
  x-averaged `f` on the beam's flank falls from 0.096 to 0.0048 by t = 120. It
  rises and falls with the swing, from 0.0013 to 0.026 over t ∈ [80, 160], so
  the test asks only for a factor of 3.
* **`f` stays a distribution function.** Its minimum is 3.6e-23, its mass holds
  to 3.4e-16, and its maximum stays below the initial one.

The seed is not the literature's 0.04: there the mode starts at a quarter of its
saturated amplitude, and there is no linear phase to fit. The default suite
checks `bump_on_tail_root` against a quadrature of the formula exactly as Arber
and Vann write it, which never calls `Z` and leaves 5.2e-15 (a beam twice as
narrow would leave 0.21). It cools the beam and follows the root to the
cold-beam fluid limit, which it meets to 0.026%, and it finds the band's edge at
k = 0.4824, where the phase velocity equals the bottom of the valley in `F` to
the last digit, as Penrose's criterion requires.

### The laser wakefield

Before the ponderomotive coupling was added, the wakefield example was asserted
only to run and stay bounded, because that was all it could support: without
the coupling the laser never entered the longitudinal push, and the wake it drew
was the slab edges relaxing. Its numbers came out bit-identical whether the
transverse current was right or wrong by thirty-two orders of magnitude. The
coupling, the force `−∂Φ/∂x` with `Φ = (pʸ² + pᶻ²)/2` in the momentum advection,
is now there, and the wakefield rows of the table under
[What the tests assert](@ref) measure the wake it produces against linear
wakefield theory rather than against a bound (`wakefield`, `linear_wake`,
`vg_pulse` and `wake_wavelength` in the harness; testset "Laser wakefield: the
laser drives a plasma wave").

The transverse momentum needed no repair, only a name. Nothing depends on `y` or
`z`, so `p⊥ + A⊥` is conserved along a trajectory and the plasma is at rest ahead
of the pulse; with `E⊥ = −∂A⊥/∂t` the accumulation `p⊥ += E⊥Δt` *is* the
canonical `p⊥ = −A⊥`, exact, not a force integral missing `v×B`. Three things
had to be fixed before the wake could be measured: `laser_amplitude` scaled
nothing, the pulse was injected by its tail, its peak due at t = −30.8, before
the run began, and `laser_duration = 5·2π` was `k_pσ = 5.7`,
which suppresses a wake by seven orders of magnitude. The defaults now are
`laser_duration = 2π`, for `k_pσ = 1.13`; `laser_amplitude = 0.3`, so that the
`1/γ` the non-relativistic `Φ` drops is worth 3.2%; `plasma_temperature = 0.01`,
to keep the edge sheaths off the measurement; and `total_time = 2π·22`, for
enough wake to fit through.

`linear_wake` integrates the driven oscillator on the `Φ` the run recorded rather
than on an idealised envelope, so the comparison assumes nothing about the shape
the laser arrives with, and it shares no code with the wake: `Φ` comes from the
FDTD side, `ex` from the Vlasov push and the Poisson solve. Its `3T` term is
load-bearing; at ten cells per wavelength, dropping the thermal correction took
the pointwise agreement from 8.9% to 26%. `Δε/ε`, 13.7%, is supposed to rise:
`ε` is the longitudinal kinetic and electrostatic energy, so a wake trading one
for the other leaves it alone, but the transverse motion and field are outside
it, and a laser doing work on the plasma pushes energy across that boundary.

The study runs at twenty cells per laser wavelength, and that is a physics
choice rather than a taste. At the ten it used before, the Yee dispersion
relation puts the carrier's group velocity 4.5% below `√(1−n)`, and the pulse,
averaged over its spectrum, 5.2% below, 6.6% as measured: the pulse arrives late
and the wake keeps station with it, so the wake's phase velocity is wrong by as
much, `γ_φ ≈ 2.2` where the physics gives 3.0. The pulse's own
bandwidth is worth only half a percent: the continuum group velocity averaged
over its spectrum is 0.9438 against 0.9487. `vg_pulse`, the group velocity of
the *discrete* relation averaged over the same spectrum, predicts 0.8993 at ten
cells against 0.8858 measured, and 0.9329 at twenty against 0.9261. None of the
assertions about the wake's *frequency* could see that, since the plasma sets
it, and the comparison against `linear_wake` cannot either, since that reference
is driven by the `Φ` of the run it is checking. It took a closed form for the
driver and a second resolution. The wake at twenty cells is 39% larger than at
ten, the most resolution-sensitive number in the study.

The `a₀²` law is tested where the signal dominates: continued to `a₀ = 0.075`
the ratio reads 3.80 and 3.33 at the two resolutions, because the unlit
control's 3.4e-4 is 43% of that wake. Two approximations remain, and the
defaults are chosen to stay inside them rather than to be impressive: the
ponderomotive potential is the non-relativistic one and the transverse current
is taken through momentum rather than velocity, both corrections of relative
order `p⊥²`, 3.2% at `a₀ = 0.3`. `test/test_em_plasma.jl` holds the transverse
half alone to the plasma's discrete dispersion relation, to 3.2e-4 over two
densities, three Courant numbers and four wavenumbers, with its cutoff and the
sign of its current. See the docstring on `wakefield`.

### The scheme comparison

The comparison of the advection schemes is made on the physics rather than on a
shifted sine, and is advisory rather than asserted
(`verification/scheme-comparison.jl`, which exits 0; the timing column is wall
clock on whatever machine runs it).

Ranked by error in the Landau damping rate at k = 0.5, `LaxWendroff` (0.67%) and
cubic `SemiLagrangian` (1.09%) lead the schemes that dissipate, ahead of PFC
(1.33%) and `Godunov` with `VanLeer` (1.45%), and upwind trails at 48.8%, with
`Godunov(PiecewiseConstant())` and linear `SemiLagrangian`, identical to it: its
own dissipation is two orders of magnitude larger than the damping it is
measuring. A 1% perturbation cannot rank the schemes on positivity, since every
one returns the same `min f = 1.3e-4`, the Maxwellian's tail at `v = ±4`. At 50%
amplitude the two that lead are exactly the two that drive `f` negative,
`LaxWendroff` to −0.093 and cubic `SemiLagrangian` to −0.105 against a peak of
0.6, over one run to t = 40 on the round trip's 64 × 121 grid: Godunov's
theorem arriving in the physics, and the reason the solvers default to PFC.
`Superbee` heads the table, at 0.21% in the rate and 0.14% in the frequency,
from below: its compression anti-diffuses, it shows the only growing
L² in the table, +8.9e-5, and its lead is a cancellation. Refined at a fixed
Courant number from Nx = 32 to 256:

| scheme | γ error | L² change |
|---|---|---|
| `Superbee` | −2.58%, −0.21%, −0.67%, +0.13% | rises at every level, +5.0e-4 to +7.3e-6 |
| `VanLeer` | +7.73%, +1.45%, +0.53%, +0.44% | falls |
| PFC | +6.47% to +0.47%, monotone | falls |

At 50% amplitude `Superbee` stays positive, at 3.05e-9.

On the kernels alone, `test/test_comparison.jl` asserts the L² error after one
traversal at N = 512:

| scheme | sine | gaussian | square |
|---|---|---|---|
| cubic `SemiLagrangian` | 2.46e-8 | 3.55e-6 | 2.01e-2 |
| PFC | 2.30e-7 | 3.30e-5 | 2.65e-2 |
| quadratic `SemiLagrangian` | 5.02e-6 | 1.37e-4 | 2.72e-2 |
| `LaxWendroff` | 4.68e-5 | 1.28e-3 | 4.53e-2 |
| `Godunov` with `VanLeer` | 8.55e-5 | 9.50e-4 | 2.77e-2 |
| `Godunov` with `Superbee` | 1.74e-4 | 1.98e-3 | 1.58e-2 |
| `Upwind` | 8.08e-3 | 4.11e-2 | 6.32e-2 |

The ranking follows the schemes' order on the sine, over five decades, and
collapses on the discontinuity to a factor of four between best and worst. The
two second-order schemes change places there: `LaxWendroff` is linear and so, by
Godunov's theorem, cannot be monotone; it rings, and falls behind the limited
`Godunov`, whose limiter buys monotonicity with the extrema it clips on the
sine. No ordering of these survives a change of problem class. `Upwind`,
`Godunov(PiecewiseConstant())` and linear `SemiLagrangian` agree to 2.9e-14
after 1280 steps, as the amplification factor they share requires
(`test/test_amplification.jl`, where six schemes are held mode by mode to the
symbols of their update formulas, to 1.8e-14, and the dissipation and phase
error per wavenumber follow from the same symbols: at λ = 4Δx the phase error is
2.4% for PFC against 29% for `LaxWendroff`).

`Superbee`, `r -> max(0, min(2r, 1), min(r, 2))` [Roe, Annu. Rev. Fluid Mech.
18, 337 (1986)], is the upper edge of Sweby's region, the most compressive
limiter that keeps `Godunov(PiecewiseLinear())` second order and total-variation
diminishing to `|c| = 1`. It is TVD from `c = 0.5` to 1 in both directions, the
pulse's total-variation ratio between 0.99999999 and 1. It is second order in
L¹ but gets there late: the local slopes from N = 32 to 2048 are 1.19, 1.81,
1.91, 1.96, 1.98 and 1.99, and `test/test_convergence.jl` holds it to 2 over 128
to 1024, where it measures 1.95. Compression has a price the TVD property does
not show: under `Superbee` the L² norm of a smooth perturbation *grows*, by
1.3e-3 on the sine and 1.2e-2 on the gaussian after a traversal at N = 128, where
every other scheme's shrinks, and `test/test_invariants.jl` asserts the sign.

`benchmark/workprecision.jl` (advisory) puts error against the cost of reaching
it and prints the frontier per problem class, the schemes no other beats on both
axes. On smooth data that is `Upwind`, `LaxWendroff`, PFC and cubic
`SemiLagrangian`, from 8.1e-3 down to 2.5e-8. On the gaussian and the pulse the
limited `Godunov` joins it: with `VanLeer` it reaches the pulse's 2.77e-2, a
little behind PFC's 2.65e-2, at under a third of PFC's cost, and with `Superbee`
1.58e-2. `Upwind` and `Godunov(PiecewiseConstant())` are the same scheme, so the
gap between their timings is the report's own error bar.

Unit conventions are in [the normalization notes](normalization.md).
