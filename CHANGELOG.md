# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project has not been released; entries below describe work on `master`.

## [0.2.0] - unreleased

### Added

- **`Superbee`**, `r -> max(0, min(2r, 1), min(r, 2))` [Roe, Annu. Rev. Fluid
  Mech. 18, 337 (1986)]: the upper edge of Sweby's region, and so the most
  compressive limiter that keeps `Godunov(PiecewiseLinear())` second order and
  total-variation diminishing to `|c| = 1`. `VanLeer` used to lead the square
  pulse in `test_comparison.jl` only because the flux lacked its `(1 − |c|)`
  factor (see Fixed). That made its limiter `VanLeer/(1 − c)`, past this edge,
  which also made the scheme first order and unstable past `c = 1/2`. Superbee
  sharpens a jump from on the edge rather than outside it.

  Measured with `c = 0.4` and N = 512 unless noted:

  * **The pulse.** The L² error after one traversal is 1.58e-2, the lowest in
    `test_comparison.jl`; the cubic spline leaves 2.01e-2 and `VanLeer`
    2.77e-2. In `benchmark/workprecision.jl` Superbee joins the pulse's
    frontier at 7.2 ms and pushes the cubic spline, at 49 ms for 2.01e-2, off
    it.
  * **Smooth data.** 1.74e-4 on the sine and 1.98e-3 on the gaussian, twice
    `VanLeer`'s 8.55e-5 and 9.50e-4. Against the cubic spline it goes from 7070
    times its error on the sine to 0.786 times it on the pulse. It is second
    order in L¹, but gets there late: the local slopes from N = 32 to 2048 are
    1.19, 1.81, 1.91, 1.96, 1.98 and 1.99, so a fit over 32 to 256 reads 1.65.
    `test_convergence.jl` holds it to 2 over 128 to 1024, where it measures
    1.95. L² settles near 1.69 and L∞ near 1.3.
  * **Stability.** It is TVD from `c = 0.5` to 1 in both directions: the
    pulse's total-variation ratio is 0.99999999 to 1, `f` stays inside
    `[0, 1]`, and a sine touching zero stays non-negative through 1000 steps.
    The mirror is exact, the step at `c = ±1` is one ulp from a shift, and it
    allocates nothing.
  * **Anti-diffusion.** Compression has a price the TVD property does not
    show. Under Superbee the L² norm of a smooth perturbation *grows*: +1.3e-3
    on the sine and +1.2e-2 on the gaussian after a traversal at N = 128, where
    every other scheme's shrinks. `test_invariants.jl` asserts the sign. In
    `verification/scheme-comparison.jl` this puts Superbee at the head of the
    Landau table, 0.39% off in the rate and 0.14% in the frequency. It gets
    there from below, with the only growing L² in the table, +8.9e-5. Refined
    at a fixed Courant number from Nx = 32 to 256:

    | scheme | γ error | L² change |
    |---|---|---|
    | Superbee | −2.41%, −0.39%, −0.68%, +0.13% | rises at every level, +5.0e-4 to +7.3e-6 |
    | `VanLeer` | +7.76%, +1.43%, +0.52%, +0.43% (monotone) | falls |
    | `PFC` | +6.41% to +0.47% (monotone) | falls |

    Its lead is a cancellation. At 50% amplitude it stays positive (3.08e-9).
  * **Cost.** 18.5 ns per cell per step on smooth data and 10.7 on the pulse,
    against `VanLeer`'s 17.2 and 9.1 in the same run.

  It joins `uniform_schemes`, so the golden, symmetry, contract and
  type-stability suites cover it. `test/data/golden.txt` gains one line; the
  other nine are unchanged bit for bit. It also joins:

  * the allocation gate;
  * the flux-limiter tests, which check Superbee piece by piece and `VanLeer`
    below it everywhere;
  * `test_1d_advection`, where it gives the donor-cell answer on the pulse;
  * the benchmark suite, which times it but does not judge it until the next
    `--rebaseline` stores a baseline;
  * both reports.

- **A nonlinear equilibrium, run to see whether it stays put**
  (`test/test_verification.jl`, `bgk_equilibrium` in the harness,
  `verification/bgk-equilibrium.jl`). Every other Vlasov–Poisson run in the
  suite starts off equilibrium and is judged on how it moves. This one starts on
  one — a function of the particle energy in the potential `−ψ cos kx`, on an ion
  background built to hold it — and is judged on how little it moves. It is also
  the first run whose physics is carried by trapped particles: at ψ = 0.5 two
  thirds of them are.

  The Maxwell–Boltzmann equilibrium, analytic across the separatrix, holds f to
  1.52e-3 of its peak and the field to 3.85e-3 through t = 50. That is not an
  oscillation about the equilibrium but a steady drift, and it is the scheme's
  dissipation: L² falls by 1.13e-3 and the entropy rises by 1.27e-4, where the
  exact flow keeps both. It converges at PFC's third order, 7.6 times less in f
  and 6.1 in the field per halving, and halving Δt alone moves neither.

  Give the trapped particles a temperature of their own and `F` keeps its value
  but not its slope across the separatrix. The worst cell is then *on* the
  separatrix at both resolutions — where the smooth case's is 13 and 24 cells
  away — twelve times the smooth error, converging at first order (1.75 per
  halving). Trapped-particle structure lives or dies there, and so does the
  accuracy of anything that has trapped particles in it.

  The teeth: ions built for a potential 10% off put the field 49% and 59% off at
  once and f drifts 24 and 27 times further; ions built on the documentation's
  old Poisson sign hold the reversed field (below, under Fixed).

  `vlasov_poisson` takes `nᵢ`, an ion density profile, and `renormalize`, which
  passing `nᵢ` turns off: the rescaling of `f` to the ions' charge is a
  trapezoid, exact only for proportional profiles, and on this matched pair it
  is 1 − 1.9e-3 and doubles the equilibrium's drift, inside its tolerances.

- **The plasma echo, with the field off and with it on** (`test/echo.jl`,
  `test/VlasovSolver/test_echo.jl`, `test/test_verification.jl`,
  `verification/plasma-echo.jl`). Phase mixing is reversible, and nothing in
  the suite checked that the solver keeps it so: free streaming measures the
  decay of a density mode, never whether the information it decayed into
  survives. A seed at k₁ = 1, a velocity kick at k₂ = 3/2 when t = 5, and around
  t = 15 a mode nobody seeded, at k₃ = 1/2.

  With the field off the echo has a closed form, exact in both amplitudes:
  `−iα·J₁(k₃ε(t − τ))·exp(−(k₃t − k₂τ)²/2)`. The Bessel function is not a
  refinement — its small-argument limit is 14.6% high at these amplitudes — and
  the sign is asserted separately, since the direction of the kick changes it
  and nothing about the peak. PFC at 128 × 241 follows the closed form to
  1.25e-3 of the peak, pointwise, out of a seeded mode the kick found at 1.84e-6:
  2.4e4 times smaller than the echo it turns into. The loss converges at orders
  2.89 and 2.78 as Δv halves, and it ranks the schemes: cubic SemiLagrangian
  1.07e-3, PFC 1.03e-2, LaxWendroff 4.48e-2, and Upwind 0.449, which returns 55%
  of the echo and loses most of the rest in the x-sweep rather than in the kick.
  Fast suite.

  With the field on, `echo_second_order` composes three linear responses of the
  Maxwellian: the seed's filament screened at k₁, the kick screened at k₂
  together with the lifetime of the field it induces, and the echo's own density
  polarising the plasma at k₃, solved as a Volterra equation. Its pieces are
  checked before any run is held to it — the resonant dielectric function is
  `dielectric` to 1.2e-16, the filament's closed form agrees with a time-domain
  solve to 4.9e-6, and with the field off the theory is the closed form's
  small-ε limit to 1.1e-14. The self-consistent run (128 × 241, α = ε = 0.01,
  τ = 10) follows it to 4.8e-3 of the peak, pointwise and ringing included, with
  the peak 0.47% low on the same sample, t = 29.05. The field-off form would put
  that peak 2.1 times higher and 1.14 later. The residual is the x-sweep's: 1.6e-2
  at Nx = 64, where halving Δv or Δt instead leaves it at 1.5e-2 and 1.7e-2.
  Extended suite.

- **Time reversal of the self-consistent system** (`test/test_verification.jl`).
  `test_strang_splitting.jl` measures reversibility on a rigid rotation with the
  field switched off; this is the same statement for Vlasov–Poisson, which is
  invariant under `(t, v) → (-t, -v)` with `E` unchanged. Run forward, flip the
  velocity axis, run forward again, flip back: an exact scheme returns the
  initial state, and what a real one loses is its own dissipation.

  Three measurements, and the second is the interesting one:

  * at α = 0.05 the round-trip error is 2.82e-3, 6.48e-4 and 1.58e-4 as the grid
    halves — ratios 4.4 and 4.1, second order, so the irreversibility is the
    scheme's and it converges away. (It was 1.97e-3, 3.44e-4 and 4.85e-5, ratios
    5.7 and 7.1 toward `PFC`'s third order, while the harness bounded the limiter
    at 1; bounded at the initial maximum, as it now is, the limiter clips the
    peak cell and that costs the order — see Fixed.);
  * at α = 0.5 it is 0.164 and stays 0.156 when the grid is halved, a factor of
    1.05. By then the flow has folded the distribution into filaments finer than
    Δv, and the information needed to run the film backwards is not on the mesh
    to be refined. Every nonlinear run in this suite is in that regime by the
    time it is interesting;
  * ranked by round trip at α = 0.5: LaxWendroff 0.085, cubic SemiLagrangian
    0.106, PFC 0.164, Upwind 0.332 — and the first two get there by driving `f`
    to −0.094 and −0.058 against a peak of 0.6, where PFC stays at +5e-10 and
    Upwind pays a fifth of its L² norm instead. Reversibility and positivity are
    the same trade-off from two sides, arrived at here through a quantity that
    has nothing to do with the damping rate `verification/scheme-comparison.jl`
    ranks them by.

  `vlasov_poisson` now returns `f` as well. It costs nothing — the array exists
  either way — and it is the only way to ask a question about the distribution
  rather than about a moment of it.

- **Strong Landau damping is measured against the literature instead of
  eyeballed.** The `α = 0.5` case has been in the notebook since 2021 with the
  text "the quantitative coincidence is almost perfect" against Fig. 6(a) of
  Filbet, Sonnendrücker and Bertrand (2001) — a comparison nothing could check.
  It is the standard nonlinear benchmark and is quoted by two numbers:
  Cheng and Knorr (1976) give γ₁ = −0.281 and γ₂ = 0.084, later work −0.292 with
  0.0815 and −0.2918 with 0.08584.

  Measured on 128 × 241 with Δt = 0.025: γ₁ = 0.2863 over the four maxima of the
  decay proper, γ₂ = 0.0789 over the eight of the regrowth. Both windows are
  conventions and the testset says so with numbers: γ₁ reads 0.3786 over three
  maxima and 0.2281 over five, because the envelope steepens and then flattens
  into the trapping plateau, so a straight line through it depends on how much
  of the curve is inside the window.

  Refinement is what makes the agreement more than a coincidence of one grid:
  γ₂ goes 0.0716 → 0.0789 → 0.0814 as Δx, Δv and Δt halve, toward the published
  value rather than away from it.

  The non-uniform velocity grid — the notebook's, and the only path this suite
  has toward an adaptive mesh — reproduces both rates to 0.04% and 0.6%. It had
  been asserted only through an energy drift on a *linear* run, where the
  distribution never approaches the sharp gradients its limiter exists for.

- **Galilean invariance of the Landau run** (`test/test_verification.jl`). Every
  Vlasov–Poisson case in this suite starts from a distribution symmetric in `v`,
  so the mean velocity is zero throughout and the drifting half of the solver
  has never been exercised — the same blind spot `test_damping_1v.jl` found in
  `BGK`, where the mean-velocity computation had never run on data with a mean
  velocity.

  A boost by `u` carries `f(x,v,t) → f(x-ut, v-u, t)` and `E(x,t) → E(x-ut,t)`,
  so the mode's history should be the same complex function times `exp(-ikut)`.
  Nothing in the discretisation is Galilean invariant — the grid does not move
  and the boosted Maxwellian sits on it asymmetrically — so the agreement is a
  measurement: over the eleven maxima in the fitting window, amplitudes within
  1.7e-3 and phases within 1.8e-3 rad, with fitted rates 0.136% apart and
  frequencies identical to the estimator's resolution. The phase is compared at
  each sample's own time, `t[k] + Δt/2` — the field is recorded mid-step — and
  removing the Doppler factor at `t[k]` instead reads 8.0e-3, of which 6.25e-3 is
  that half-step and not the grid. Without the Doppler
  correction the phase differs by up to 2.88 rad, which is what says the
  correction is doing work.

  Compared at the maxima rather than pointwise: both `|E_k|` and `ε_e` pass
  through deep nulls, and the first version of this test reported a 900%
  discrepancy that was entirely two nulls landing a time step apart.

### Fixed

- **`two_stream(0.4)` ran on 95 cells rather than 96.** `two_stream` built its
  grid as `collect(Δx:Δx:L)`, and a floating-point range works its length out
  from its endpoints: at `a = 0.4`, `Δx + 95Δx` rounds past `L` and the range
  stops a point short. The box was then `95Δx = 46.63` against `L = 47.12`, its
  fundamental `k·96/95` — `a = 0.4042`, where the warm root is 0.30536 rather
  than 0.30362 — and the seeded `cos kx`, a 96-point period on 95 points, left
  0.76% of its amplitude outside that fundamental. The −0.39% the test reported
  against the root at 0.4 was −0.95% against the root of the wavenumber the run
  actually had.

  The grid is now `range(Δx; step = Δx, length = Nx)`, `Nx` points by
  construction. On 96 cells the rate is 0.30173: −0.62% from the warm root and
  −2.10% from the cold, fitted over t ∈ [15.50, 21.60], inside the 3% tolerance
  and still below the peak at a = 0.6. What `growth_rate` and `two_stream` quote
  at a = 0.4 moved with it: fixed time windows give 10.65%, 6.43%, 3.99% and
  0.19% from the cold root (were 9.75%, 5.63%, 4.35% and 0.55%); the spread over
  start points is 41.4% over one beat period and 9.4% over two (were 39.6% and
  9.0%); `cos ω₊t` and `sin ω₊t` in the band fit move it from −2.10% to −3.24%,
  1.1 points where the note said at most one; continued past its fit, the
  velocity sweep first passes its Courant limit at t = 23.25 (was 23.20), and it
  is at 0.528 when `ε_e` first reaches 5.0 (was 0.524), which takes the range
  the Courant-bound entry below quotes for the fits to 0.53 to 0.73; and the
  growth plot's single-γ line peels 6.0× below the curve (was
  6.1×).

  Every other periodic grid in the extended suite and in
  `verification/scheme-comparison.jl` was built the same way and had its length
  right. They are built by length now too, and are the same grids to the bit —
  checked on Julia 1.10 and 1.13 — so nothing else moves. `[j*Δx for j = 1:Nx]`,
  the obvious alternative, rounds `jΔx` once where the range rounds
  `Δx + (j − 1)Δx` twice, and would have moved points of every grid here but the
  two stable two-stream cases by an ulp. The colon form is a trap rather than a
  one-off: at the two-stream construction it comes up short for 14 of the 156
  values of `a` in `0.05:0.01:1.6`.

  One quote was wrong on its own account. `growth_rate`'s note credited fitting
  "over an integer number of beat periods instead of an arbitrary window" with
  the drop in that spread; the measurement behind it compared one beat period
  with two, and the note now says so.

- **`Godunov(PiecewiseLinear())` is second order, and total-variation
  diminishing to `|c| = 1`.** Its flux was the reconstruction's value at the
  interface, `|c|(fᵢ₋₁ + φ(r)(fᵢ − fᵢ₋₁)/2)`, which makes the update forward
  Euler on a limited slope. It is now the reconstruction averaged over the strip
  that crosses the interface in one step. That strip's midpoint lies `|c|Δx/2`
  upwind of the interface, `(1 − |c|)Δx/2` downwind of the cell's centre, so the
  slope term gains a factor `(1 − |c|)`, and the scheme becomes what its name
  says: Godunov's
  reconstruct–evolve–average with a linear reconstruction, which is Sweby's
  flux-limited Lax–Wendroff. The old flux is Sweby's as well, with the limiter
  `φ/(1 − c)` in place of `φ` (the two agree to 6.7e-16 after 50 steps), and
  that one substitution accounts for all of its behaviour:

  * **It was stable only to `|c| ≤ 1/2`**, the limit to which `φ/(1 − c)`
    satisfies Harten's condition. With `VanLeer` at N = 128, 200 steps
    multiplied a square pulse's total variation by 2.05 at `c = 0.6`, 39 at
    0.7, 4.8e5 at 0.8 and 1.0e7 at 0.9. `0.5(1 + sin)` went negative in one
    step from `c = 0.58`, to −9.0e-6. `1 + 0.5 sin` was 84 from the exact
    answer after 300 steps at 0.7, and 893 after 100 at 0.9. Now the pulse's
    total-variation ratio is 0.999992 to 1 at every `c` from 0.5 to 1 in both
    directions, `f` stays inside `[0, 1]`, and the sine stays non-negative
    through 1000 steps (`test_invariants.jl`). In the same two runs
    `1 + 0.5 sin` is 3.5e-3 and 1.1e-3 from the exact answer.
  * **It was first order.** `φ/(1 − c)` is `1/(1 − c)` at `r = 1`, where second
    order needs 1, so it steepened every smooth slope. At `c = 0.4` the L¹, L²
    and L∞ orders were 1.00, 1.01 and 0.82; the last had been quoted as 0.69.
    They are now 2.09, 1.77 and 1.42, and `test_convergence.jl` holds the L¹
    order to 2 where it held it to 1.
  * **Without a limiter it was unstable at every `c`.** `NoLimiter`'s 1 made
    the flux centred, `|g| = √(1 + c²sin²θ)`, and it reached 1.4e25 times its
    amplitude in four traversals at `c = 0.4`. It is now `LaxWendroff`: 2.2e-16
    from it after one step and 1.4e-14 after 1280, and `test_amplification.jl`
    checks it against Lax–Wendroff's symbol mode by mode.

  At `|c| = 1` the factor vanishes and the step is an exact shift. It is now one
  ulp from `circshift`; it was 2.4e-3 off without a limiter and 2.7e-3 with
  `VanLeer`, and `test_contracts.jl` no longer excuses it. The `|c| ≤ 1` that
  `advect!` enforces is now this scheme's real limit.

  The same substitution is what made the old scheme good at a jump. From
  `c = 1/3` up, `VanLeer/(1 − c)` lies above Superbee, the most compressive
  limiter inside Sweby's region, at every `r`, and it steepened the pulse just
  as it steepened the sine. With `VanLeer` at N = 512 and `c = 0.4`, the L²
  error after one traversal (`test_comparison.jl`) moves as follows:

  * sine: 7.03e-3 → 8.55e-5;
  * gaussian: 1.61e-2 → 9.50e-4;
  * square pulse: 8.45e-3 → 2.77e-2, a rise. The scheme led that column and now
    sits fourth of six.

  Superbee itself, which keeps the scheme second order and total-variation
  diminishing to `|c| = 1`, measures 1.58e-2 on the pulse; it is now a limiter
  of its own (see Added). In the Landau damping comparison
  (`verification/scheme-comparison.jl`) the rate error goes from 3.99% to 1.43%
  and the frequency error from 16.0% to 0.54%. In the work–precision report the
  scheme is now dominated on every profile, by `PFC`.

  The factor costs one subtraction and one multiplication per flux. The
  `VanLeer` step, 17 ns per cell on the sine and 9 on the pulse, measured 2% to
  9% slower with it. Those are medians of three alternations at N = 512 and
  10000, and the step was slower in every configuration. `LaxWendroff` and
  `Upwind`, whose code is untouched, moved by −11% to +10% over the same runs.
  Without a limiter the step costs 0.4 ns per cell and did not move measurably.

  `test/data/golden.txt` is regenerated for `Godunov_linear` and
  `Godunov_linear_VanLeer`. The other seven lines are unchanged bit for bit.

  Found on the way: `test_amplification.jl` quoted its three-point schemes as
  agreeing with their symbols to 8.9e-15. They agree to 1.6e-14 to 1.8e-14, on
  Julia 1.10.12 and 1.13.0 alike, and the quote is corrected.

- **`advect!` refuses a step past the scheme's Courant limit.** `_validate`
  checked aliasing, the lengths and the workspace, and not the one bound every
  explicit scheme here has. `Upwind`, `LaxWendroff`, `Godunov` and `PFC` are
  unstable past `|c| = 1`, and they fail plausibly, because what grows is
  round-off at the grid scale. One step on `1 + 0.5 sin` at N = 128, against the
  exact shift, with the worst per-step growth over every mode of the grid and
  the error a hundred steps later:

  | scheme | c = 0.9 | 1.0 | 1.2 | 2.0 | growth at 1.2 | 100 steps at 1.2 |
  |---|---|---|---|---|---|---|
  | `Upwind` | 5.4e-5 | 0 | 1.4e-4 | 1.2e-3 | 1.40 | 2.7e-2, and 2.0e27 at 300 |
  | `Godunov(PiecewiseConstant())` | 5.4e-5 | 2.2e-16 | 1.4e-4 | 1.2e-3 | 1.40 | 6.8e-2, and 5.2e27 at 300 |
  | `LaxWendroff` | 1.7e-6 | 2.2e-16 | 5.2e-6 | 5.9e-5 | 1.88 | 9.1e10 |
  | `Godunov(PiecewiseLinear(), VanLeer())` | 1.4e-5 | 2.2e-16 | 3.6e-5 | 3.0e-4 | (nonlinear) | 2.2e10 |
  | `PFC` | 2.3e-8 | 2.2e-16 | 5.1e-8 | 8.9e-16 | 1.18 | 5.1e-6, and 4.9e15 at 300 |
  | `SemiLagrangian` cubic | 9.8e-10 | 2.2e-16 | 3.1e-9 | 2.2e-16 | 1.00 | 3.1e-7 |

  The `Godunov(PiecewiseLinear(), VanLeer())` row is for the flux as corrected
  in the entry above, measured past `|c| = 1` on a copy of the kernel that
  matches it bit for bit where `advect!` accepts the step. Before the
  correction the row read 5.5e-4, 6.8e-4, 9.6e-4, 2.6e-3 and 1.6e13.

  `PFC` loses first the property it exists for: one step at 1.2 takes
  `0.5(1 + sin)`, which touches zero, to −2.9e-5, and to 1 + 2.9e-5 against an
  `fmax` of 1. Its exact answer at `c = 2` is a coincidence of the unlimited
  reconstruction, which interpolates the primitive at the stencil's nodes and so
  reproduces a whole-cell shift; at 1.5 it goes negative in one step and reaches
  5.6e276 in a thousand. A free-streaming run whose fastest rows sat at
  `c = 1.22` (`vmax = 6`, `Δt = 0.02`, `Δx = 4π/128`) returned silently from
  every such call, and the first sign was `PFC`'s own `checked` assertion — 218
  steps in, at `minimum(src) = −1.13e-10`, in a reproduction. With
  `checked = false`, or with `LaxWendroff` or `Upwind`, there would have been
  none.

  `|c| > 1` is now a `DomainError`, raised before anything is written, and so is
  `NaN`. `c = ±1` is accepted: it is an exact one-cell shift for every bounded
  scheme. The bounded method is the default, and
  `SemiLagrangian`, whose characteristic tracing has no Courant limit, opts out,
  so a scheme added later is checked unless it says otherwise. `PFCNonUniform`
  takes a displacement and is held to its narrowest cell, stored at
  construction: every cell gives up its flux alone, and on a 1:2 grid a step of
  1.1 narrow widths — only 0.55 of the wide one — takes data with an exact zero
  to −7.6e-6. The check does not depend on `checked`. It is one comparison
  outside the loop, under 3 ns as a call on its own — 0.16% of `Upwind`'s 1.9 µs
  step at N = 10000 — and inside the step it is lost in the noise of timing it:
  with it and without, alternated eight times at N = 10000, the schemes ran
  0.3% to 7.3% faster *with* it, against a 2.2% shift in a scheme it does not
  touch. `PFC`'s `checked` pass, timed alongside, reads 8.0% where its docstring
  quotes 13%.

  Two runs in the extended suite, and the notebooks that mirror them, had been
  taking such steps; a probe on every `advect!` call of both suites and every
  notebook found no others, the largest being 0.94 in the `a = 1.0` two-stream
  run. The Landau x-sweeps sit at 0.32 to 0.41, being Strang half-steps; it is
  the velocity sweeps, full steps driven by the field, that crossed:

  * the two-stream runs keep growing after their fits, and the velocity sweep
    reached 217 cells a step before the `a = 0.6` run went to `NaN` at
    t = 24.15, which is what bounded `tmax` from above;
  * the strong-damping run on the notebook's non-uniform velocity grid starts
    with the field at 1.0017 and `Δt` equal to the narrow cells' width, so its
    first step asks for 1.0017 of them.

  `line_advector` now splits a displacement wider than the narrowest cell into
  the fewest sub-steps that fit, and `verification/landau-damping-1d1v.jl` does
  the same. No fitted number moves: every two-stream fit ends at a velocity
  Courant number of 0.53 to 0.73, before the first split, so the rates are
  bit-identical, and the strong-damping run's four split calls move γ₁ and γ₂ in
  the sixth digit. The two-stream runs now saturate instead of diverging — `ε_e`
  peaks at 106 for `a = 0.6` — and are finite to t = 40, so `tmax` has no upper
  bound left. `test_invariants.jl` marched `Upwind` at c = 1.05 for 2000 steps
  to show the instability; it asserts the refusal instead.

  Measured on the way: `Godunov(PiecewiseLinear())` was stable only to
  `|c| ≤ 1/2`, because its flux lacked a `(1 − |c|)` factor. That is fixed in
  its own entry above, and `|c| ≤ 1` is now its limit as well.
  And the velocity Courant number `growth_rate` quoted for the two-stream fits,
  0.46 to 0.65, was the single-mode estimate `√(2ε_e/L)·Δt/Δv`; the largest over
  `x` is 0.53 to 0.73.

- **A full test run evaluates each shared test file once.** `runtests.jl`
  includes every test file into `Main`, and the files several of them share were
  evaluated once per file that included them: `scheme_cases.jl` ten times,
  `dispersion.jl` six, `echo.jl` three and `verification_harness.jl` twice. Each
  repeat redefined the file's methods, and `Pkg.test` runs with
  `--warn-overwrite=yes`: 150 warnings, in the default run and the extended one
  alike, among which one that mattered would not have been read. Each include of
  a shared file in `test/` is now guarded on a function that file defines —
  `@isdefined(march!) || include(…)`, the idiom Julia's own test suite uses for
  its helpers — so every test file still runs on its own and a full run prints
  none of the 150. Guarded on a name rather than on the path, because an
  `include_once` helper would itself have to be included, and guarded, by every
  file that used it. The name is one a session would not have for itself:
  `dispersion.jl` is guarded on `landau_root` rather than `Z`, whose skip made
  `test_dispersion.jl` test a session's own `Z(x) = x` — 12 failures — where the
  file would have replaced it. The scripts in `verification/` and `benchmark/`
  keep their plain `include`: each is the first thing in its process to include
  the file. The harness's own includes are guarded as well, so a script that
  includes it now evaluates `dispersion.jl` once where it did twice.
  `test_allocations.jl`, which included `scheme_cases.jl` but keeps its own list
  — the canonical one less the two spline schemes, whose prefilter allocates and
  is bounded apart — no longer includes it. Pass and broken counts are
  unchanged: 1499 and 2 by default, 1632 and 2 extended.

- **`test_free_streaming.jl` replaced the harness's `mode_amplitude` for every
  file that ran after it.** Its own `mode_amplitude(n, x, k)`, the real `cos kx`
  projection, had the signature of the harness's complex one, and a full run
  puts both files in `Main`: from there on, a mode the harness measured would
  have come back as its real part alone. Nothing measures one after it yet; the
  warning `Pkg.test` printed for it was the one of 151 that was not a repeat. The
  free-streaming helper is `cos_amplitude` now, and neither run prints a
  method-overwrite warning.

- **`PFC` in the remaining verification runs is bounded by the distribution it
  carries.** The harness took its own defaults' bounds from `f` above; the runs
  that pass `PFC` in explicitly still bounded it at 0.5 (`test_echo.jl`,
  `test_free_streaming.jl`, `verification/plasma-echo.jl`) or 1
  (the reversibility ranking in `test_verification.jl`,
  `verification/scheme-comparison.jl`). A bound taken from `f₀` by the caller
  cannot work where the driver rescales `f` before running — by 0.8% at α = 0.5,
  above the bound, into `PFC`'s own check — so `vlasov_poisson`, `ballistic_echo`
  and `free_stream` now also take a scheme as a function of the initial `f`,
  `f -> PFC(fmin = 0.0, fmax = maximum(f))`, and call it with the `f` they start
  from.

  The numbers barely move: the
  echo's PFC error goes from 1.23e-3 to 1.25e-3 of the peak (and 9.83e-3 to
  1.03e-2 at 64 × 121), its orders from 2.87 and 2.72 to 2.89 and 2.78, PFC's
  L² loss in the α = 0.5 round trip from 0.047 to 0.048, and in the scheme
  comparison the uniform `PFC` row now equals `PFCNonUniform`'s, 1.31% where it
  read 1.45%. Free streaming does not move at all. Left as they are, on purpose:
  the unit tests of the schemes, where a bound is an input of the contract under
  test rather than the maximum of an evolving distribution, and the
  rigid-rotation test, whose Gaussian's analytic peak is the 1 it is given.

- **`docs/normalization.md` gave the Poisson equation with the wrong sign.** It
  read `∂E/∂x = nᵢ − nₑ` beside `∂f/∂t + v∂f/∂x + E∂f/∂v = 0`, a pair under which
  electrons attract one another. Every Poisson call in the package and the
  verification runs solves `∂E/∂x = nₑ − nᵢ`, `E` being the force on an electron
  — minus the physical field — as the transverse section of the same page already
  had it. It was found by building on it: ions made for a nonlinear equilibrium
  under the documented sign hold exactly the reversed field, `E/E₀ = −1.000` from
  the first sample, and the run leaves the equilibrium by 16% of its peak.

- **The verification harness gave PFC an upper bound of 1 whatever `f` was.**
  `vlasov_poisson` built its default `PFCNonUniform` schemes with `fmax = 1.0`,
  a number nothing chose, and `PFCNonUniform` does not check its data against
  it. A trapped population at half the passing temperature peaks at 1.79, and
  the limiter took it 43.6% of its peak away from equilibrium by t = 50. The
  bounds are now those of `f` on entry, 0 and its maximum, which by Liouville
  the exact solution keeps: that run holds to 0.99%. The same goes for the
  other two places the harness chose a bound — the echo's kick (1) and
  `two_stream` (3, set once to clear the 1) — and for the Landau and
  plasma-oscillation notebooks, which the harness mirrors. No run exceeds its
  bound: a new `fmax` history tops out exactly at it or below, to the bit.

  Every earlier run sat under the old bound, where the limiter never engaged.
  At the maximum it now does, and the price is measured: the α = 0.05
  reversibility round trip converges at second order, 2.82e-3 → 6.48e-4 →
  1.58e-4, where it went 1.97e-3 → 3.44e-4 → 4.85e-5 at third, and its
  threshold is now 3.5× per halving rather than 4×. The rest moved little —
  the three Landau rates by 0.2% at most (0.5% on the coarsest level of the
  refinement ladder), ω at k = 0.5 from 0.14% to 0.26% off the root, the
  two-stream rates in the fifth decimal place — and every quoted number is updated
  to the new runs. Six that were already stale are corrected with them: the
  trapezoid's mass drift in the invariants note (1.6e-4 quoted, 2.4e-4 measured
  every step, under either bound); the recurrence window ratio (116 quoted, 119
  at the old bound); the stable two-stream cases' decay (3.36e-5 and 7.02e-6
  quoted, 2.54e-6 and 1.2e-9 measured) and the unstable case they are set
  against; the `a = 1.0` row of the warm-root table (≈0.088, from before that
  case ran to `tmax = 80`, against 0.09516); and the recurrence table's floor
  and `two_stream`'s "peak ε_e of 61", which no definition reproduces and which
  are replaced by measured quantities.

- **The "linear" Landau case at k = 0.3 was not linear, and the agreement it
  reported was two errors cancelling.** At the 1% perturbation this suite used,
  the local damping rate leaves the analytic value at t ≈ 35, reaches zero at
  t ≈ 73 and goes *negative* after that: the field grows again as trapped
  particles slosh. The fit ran to t = 50, inside that, and came out 0.42%
  **below** the analytic rate — the only measurement in the suite that sat low,
  where numerical dissipation can only push high.

  The note on file said the flattening was "the estimator running out of signal,
  not the physics changing", and offered as evidence that it does not move when
  Δx, Δv and Δt are halved. That is evidence for the opposite: a numerical
  artefact moves under refinement and physics does not.

  At α = 1e-3 the same column is flat to t = 95, the window widens to [10, 90]
  with thirty maxima, and the measurement reads 0.71% *above* — as do the other
  two, which is what a dissipative scheme should do. The three cases now run at
  α = 1e-3 and are 0.71%, 1.14% and 0.94% on γ, 0.08%, 0.22% and 0.26% on ω,
  with window spreads of 0.18%, 0.05% and 0.04% (they were up to 1.45%).

  `trapping_phase` is the guard, asserted per case: `(√α/γ)(1 - exp(-γT))`, the
  bounce phase accumulated before the mode damps away, evaluated at the
  amplitude each run actually used. It is calibrated by the new trapping testset
  rather than assumed: in its own units the damping departs at about 2.8 and
  stops at 4.3 to 5.7, and the three cases sit at 0.21, 0.45 and 1.70 where the
  old k = 0.3 setup sat at 3.7. (The testset's headline invariant is the
  undamped `√α·t₀`, 7.1 to 8.0 at the arrest — a different number, not to be
  compared with these.)

- **The recurrence in the Landau notebook is the second harmonic's, not the
  seeded mode's.** The text explained the rise in `ε_e` at t ≈ 62 as the mode
  returning at "π/(kΔv)". The right expression is `2π/(kΔv)` = 125.7; what
  arrives at 62.8 is the k = 1 harmonic, generated nonlinearly, recurring at
  half the time because its wavenumber is twice as large.

  `vlasov_poisson` gained a `modes` keyword for this — complex field amplitudes
  per mode, since `ε_e` sums the box and cannot tell two modes apart. Asserted:
  the seeded mode peaks at t = 128.6 against 125.7 and comes back with 52% of
  its amplitude, the harmonic at 64.3 against 62.8 from a floor of 2.6e-17, and
  at t = 64 the harmonic is 42 times the mode that was seeded. The notebook runs
  to t = 140 now so that both recurrences are inside it, plots them separately,
  and takes its theory curve from `landau_root` — the asymptotic it plotted
  before is a rate for the energy, 2γ, which is why it looked right against
  `ε_e` while being twice γ.

- **The laser in the wakefield study was travelling six percent slow, and the
  explanation on file was wrong.** `wake_wavelength` said there was no closed
  form to predict the driver with, because a one-cycle pulse "has a bandwidth of
  order its own carrier and no single ω₀ describes it" — the measured 0.886
  against a continuum `√(1-n)` of 0.949 was put down to that.

  Measured, the bandwidth is worth half a percent: averaging the continuum group
  velocity over the pulse's own spectrum gives 0.9438 against 0.9487. The other
  six percent was the Yee grid at ten cells per wavelength. `vg_pulse` — the
  group velocity of the *discrete* relation, averaged over the same spectrum —
  predicts 0.8993 there, against 0.8858 measured; at twenty cells it predicts
  0.9329 against 0.9261.

  So `Δx` halves, to twenty cells per laser wavelength, at about four times the
  wall clock: `Δt` halves with it, so both `Nx` and `Nt` double.
  It matters because the wake keeps station with the driver: a driver six
  percent slow is a wake whose phase velocity is six percent wrong, `γ_φ ≈ 2.2`
  against 3.0, which is a different statement about trapping and dephasing.

  Nothing in the testset could see it. The wake's frequency is set by the
  plasma, not the laser's grid, and moves by 0.02% between the two resolutions;
  the comparison against `linear_wake` is blind by construction, since that
  reference is driven by the `Φ` of the run it is checking, so both sides shift
  together. What was needed was a closed form for the driver and a second
  resolution, and the testset now has both: the pulse speed is asserted against
  `vg_pulse` at each grid, and refining is asserted to move it toward the
  continuum.

  Re-measured at the new default: ω 0.23% from Bohm–Gross (was 0.36%), λ 0.37%
  from the measured driver and 1.12% from the predicted one, phase locking
  0.58%, amplitude 6.4% under linear theory with pointwise rms 6.4%, the `a₀²`
  ratio 4.13, the unlit control 36× down, `Δε/ε` 13.7%. The wake itself is 39%
  larger than the coarse grid gave — it was the most resolution-sensitive number
  in the study and nothing had been comparing it across grids.

  The `a₀²` ladder gained a rung and lost a claim: continuing to `a₀ = 0.075`
  gives ratios of 3.80 and 3.33 at the two resolutions, because the unlit
  control's 3.4e-4 is 43% of that wake. The law is tested where the signal
  dominates, and the few percent left at the top is recorded rather than
  explained.

### Added

- **The dispersion relation of an electromagnetic wave in plasma is asserted**
  (`test/test_em_plasma.jl`). The Yee update is pinned in vacuum mode by mode
  against a closed form; the plasma current that turns it into half of a
  laser-plasma model was exercised only through `wakefield`, where it enters a
  wake dominated by the ponderomotive term. Standing PEC modes in uniform
  plasma now measure

      (2/Δt)²·sin²(ωΔt/2) = (2/Δx)²·sin²(kΔx/2) + n

  over two densities, three Courant numbers and four wavenumbers, worst
  departure 3.2e-4 — which is the crossing estimator's own Δt resolution, not
  the scheme. The plasma cutoff comes with it: at n = 0.3 the longest mode the
  grid holds oscillates at 0.55018 against `√(n + k²)` = 0.55000, where the
  discrete plasma frequency itself is 0.54789 and `√n` is 0.54772.

  The sign of the current is the point of the last testset. `docs/normalization`
  says reversing it turns the oscillation into growth; run, it reaches 4.8e+134
  in four thousand steps against the 1.0e+2 the physical sign holds.

  The three lines this tests — the canonical `p⊥ += E⊥Δt`, the `-J Δt`
  convention and that sign — are now `transverse_step!` in the harness, called
  by `wakefield` rather than written inside it, so the relation asserted here is
  a statement about the code the study runs. Numbers unchanged: the extended
  suite reproduces every wakefield figure bit for bit.

  `em_omega` is the closed form. It is not the continuum `ωₚ² + k²`, and the gap
  matters: at the ten cells per wavelength the wakefield study used, the group
  velocity it implies is 4.5% below `√(1-n)`.

- **The kinetic dispersion relation is solved rather than tabulated**
  (`test/dispersion.jl`). The plasma dispersion function is written as
  `Z(ζ) = i√π·erfcx(-iζ)`, i.e. through the *entire* Faddeeva function, so one
  expression covers the whole plane: the Cauchy integral that defines `Z`
  converges only above the real axis, and every root worth finding — damped
  Landau modes, growing beam modes — lives off it. On top of that sit the
  Maxwellian susceptibility, a secant root finder, `landau_root(k)` and
  `two_stream_warm(a; vt)`.

  `test/test_dispersion.jl` holds it to identities and independent quadratures
  rather than to stored numbers: `Z` against its own Cauchy integral where that
  integral converges (3.5e-15), `Z′` against a central difference (9.3e-11,
  which is the difference's truncation and not the identity's), the
  susceptibility against a trapezoid that never calls `Z` (1.3e-14). The three
  Landau roots this suite carried as constants are now the fixture rather than
  the target, and the solver reproduces every digit they quoted.

  `SpecialFunctions` is a test-and-verification dependency for `erfcx`; the
  package itself still depends on four packages.

- **The two-stream case is measured against warm beams.** The runs have
  Maxwellian beams at `vt = 0.3` and were compared against the cold closed form,
  which at that temperature is off by up to 3.14% — half of a 6% tolerance spent
  on a known approximation, and 7.44% at the `vt = 0.6` run, where the cold error
  has grown with the temperature. Against `two_stream_warm` the same three
  measurements read −0.62%, −1.95% and +0.09%, and the tolerance is now 3%.

  `a = 1.0`, the cold stability boundary, turns from a qualitative case into the
  sharpest one in the testset. The cold form predicts exactly zero there and the
  old assertion could only say "something grows" (a factor of 4.91 over t ≤ 26);
  the warm root predicts 0.09823 and the run gives 0.09516, which is 3.12%. It
  costs a longer run — `tmax = 80` to bring `ε_e` up to the same amplitude band
  every other case uses, against a divergence at t = 86.2 — and is held to 8%
  because the beat ripple `growth_rate` documents is worst where `γ` is
  smallest: over `hi` ∈ {1, 2, 3, 5} the fit moves between −6.01% and −3.12%.

  The temperature comparison at `a = 0.6` is now absolute as well as relative:
  γ falls 4.4% between `vt = 0.3` and `vt = 0.6` where warm theory predicts
  4.4%, and each run is separately within 3% of its own root.

### Fixed

- **"Finite temperature cannot make a beam grow faster than a cold one" is
  false**, and it was load-bearing: it is why the +0.31% overshoot at
  `a = 0.8` was attributed to the beat ripple rather than to physics, and why a
  sweep in `vt` that crossed the cold value was treated as an artefact of the
  estimator. The warm root crosses *above* the cold one between `a = 0.75` and
  `a = 0.8`, and past `a = 1` the cold branch is identically zero while warm
  beams still grow — which the same testset had been measuring all along, at
  `a = 1.0`, without the two statements being put side by side. Against the warm
  root the `a = 0.8` measurement is +0.09% and there is nothing left to explain.

  The estimator's ripple is real and the note in `growth_rate` about it stands;
  what changed is that it is no longer asked to account for a discrepancy that
  belongs to the theory the measurement was being compared against.

### Added

- **The laser wakefield study drives a wake, and the test measures it against
  linear theory.** `wakefield` had **no ponderomotive coupling**: the laser
  never entered the longitudinal push, so `ex` was the slab edges relaxing, and
  the test could assert nothing beyond "runs and stays bounded". The coupling is
  the force `−∂Φ/∂x` with `Φ = (pʸ² + pᶻ²)/2` in the momentum advection, added
  alongside `e` before the sweep.

  The transverse momentum needed no repair, only a name. Nothing depends on `y`
  or `z` here, so `p⊥ + A⊥` is conserved along a trajectory and the plasma is at
  rest ahead of the pulse; with `E⊥ = −∂A⊥/∂t` the accumulation that was already
  there, `p⊥ += E⊥Δt`, *is* the canonical `p⊥ = −A⊥`. Read as a force integral
  it would be missing `v×B` and the convective term; read as the invariant it is
  exact, and the docstring now says which.

  **Three things had to be fixed before any of it could be measured.**
  `laser_amplitude` scaled nothing — the pulse shape was hard-coded to unity and
  the parameter reached only the momentum grid, so the field was bit-identical at
  `a₀` = 0.5 and 1.0. The pulse was injected by its tail: the profile was written
  in `x - t` alone, putting its maximum at `t = x_min + Δx = −30.8`, before the
  run began, and the 0.383 this study reported as its peak laser field was
  `exp(−0.96)` of the amplitude it never reached. And `laser_duration = 5·2π` was
  `k_pσ = 5.7`, which suppresses the wake by seven orders of magnitude: even with
  a coupling term, that pulse would have driven nothing.

  What is asserted now, at the study's own resolution, is that the wake
  oscillates at the Bohm–Gross frequency (0.36%), that its wavelength is the one
  a driver at the measured pulse speed imposes (1.16%), that its phase velocity
  is that pulse speed — two measurements along two different axes agreeing to
  0.75% — that its amplitude and pointwise profile match the driven-oscillator
  solution `linear_wake` integrates on the recorded drive (4.4%, and 8.9% of the
  theory's own rms), and that it is the laser's: `a₀²` to 0.45%, thirty-two
  times the unlit control, and seven and a half times larger behind the pulse
  than ahead of it.

  `linear_wake` is driven by the `Φ` the run recorded rather than by an idealised
  envelope, so the comparison assumes nothing about the shape the laser arrives
  with or whether it translates rigidly; it shares no code with the wake it is
  compared against, `Φ` being the `FDTD1D` side and `ex` the Vlasov push and the
  Poisson solve. Its `3T` term is load-bearing — dropping the thermal correction
  alone takes the pointwise agreement from 8.9% to 26% and fails the test.

  The defaults moved with the physics: `laser_duration = 2π` for `k_pσ = 1.13`,
  `laser_amplitude = 0.3` so that the `1/γ` the non-relativistic `Φ` drops is
  worth 3.2% rather than the 25% it would be at 1.0, `plasma_temperature = 0.01`
  to keep the edge sheaths off the measurement, and `total_time = 2π·22` for
  enough wake to fit through. `Δε/ε` is 8.4% against the 1.2% it was, and that is
  the point: `ε` is `∫∫f p² + ∫e²`, longitudinal kinetic *and* electrostatic, so
  a wake trading one for the other leaves it alone — but the transverse motion
  and the transverse field are outside it, and a laser doing work on the plasma
  pushes energy across that boundary. It is supposed to rise.

### Fixed

- **The two-stream growth rate appeared to overshoot the cold limit at small
  beam temperature, and does not.** The sweep quoted alongside the test was
  measured with a fixed time window of `t ∈ [10, 20]` — the exploratory
  estimator, not the amplitude band the test actually uses — and it crossed
  zero, reading +0.40% and +0.56% at `vt` = 0.2 and 0.15. Finite temperature
  cannot make a beam grow faster than a cold one, so the crossing was an
  artefact either way; under the estimator in use the sweep is monotone and
  entirely below the cold value.

  **The mechanism.** `ε_e` is not one exponential. The quadratic behind
  `γ_cold` has four roots — the growing pair `±iγ` and an oscillating pair
  `±ω₊` — and an initial perturbation excites all of them. The cross term puts
  a ripple on `ε_e` at `ω₊`, and its size *relative to* the growing mode falls
  only as `exp(−γt)`, so it survives any window one can afford. Measured at
  `a = 0.6`, the instantaneous rate oscillates with period 4.5 against the
  `2π/ω₊ = 4.626` predicted, swinging between 0.21 and 0.41 about a `γ_cold` of
  0.353. Changing `vt` changes `γ` slightly, moving the ripple's phase within a
  fixed window and dragging the fitted rate across the cold value with it.

  Fitting over two beat periods instead of one cuts the spread over start
  points from 41.4%, 14.4% and 22.3% to 9.4%, 5.6% and 4.8%. The amplitude band
  already spans 1.15 to 2.22 periods and needs no such help: adding
  `cos ω₊t`/`sin ω₊t` to the design matrix — still linear, `ω₊` being in closed
  form — moves its results by 1.1 points at most and was not kept.

  Three other explanations were measured and rejected: refining `Δv` moves the
  result by 1e-5; the driver's renormalisation leaves the effective density at
  1.0000158, worth 0.0008% on `γ`; and although the second harmonic at
  `a = 0.4` really is the more unstable of the two (`γ(0.8) = 0.311` against
  `γ(0.4) = 0.308`), it starts at `O(α²)` and gains 13% over the run against a
  head start of 1e-6.

  The test now asserts the relative statement — a colder beam grows faster at
  fixed wavenumber — rather than the tidier "every rate lies below `γ_cold`",
  which is false: `a = 0.8` comes out 0.31% above. The residual ripple biases
  either way depending on how much of a beat period the band leaves unaveraged,
  so the sign at any single wavenumber is not a property worth asserting.
  Comparing two temperatures at the same wavenumber holds the band fixed and
  leaves only the physics.

### Added

- **The two-stream instability** (`test/test_verification.jl`), the first
  *unstable* case in the suite. Everything else here is a damped or neutral
  mode, and a growth rate catches a class of error that damping cannot: a sign
  flip in the field push turns damping into growth and growth into damping, so a
  suite made only of damped cases is half-blind to it.

  The test plan deferred this twice on the grounds that a growth rate means
  either a plasma dispersion function or a hard-coded constant of unknown
  provenance. That turns out not to apply to the **cold** case, which is why it
  was worth waiting for. For two beams of density 1/2 at `±v₀` the relation
  `1 = ½/(ω − kv₀)² + ½/(ω + kv₀)²` is, with `a = kv₀` and `u = ω²`, a
  quadratic — `u± = [(2a² + 1) ± √(8a² + 1)]/2` — unstable exactly when `a < 1`,
  with `γ = √(−u₋)`. No special functions, no numerical root, nothing quoted.
  The test checks the closed form against the relation it came from (residual
  1.4e-14) and confirms it reproduces `√(3/8)` and `1/(2√2)` on its own before
  using it.

  Measured at three wavenumbers, with the beams at `vt = 0.3`: γ = 0.30173,
  0.34228 and 0.31229 against 0.30819, 0.35339 and 0.31134 — 2.10%, 3.14% and
  0.31%. **`γ(a)` is non-monotone**, peaking at `a = √(3/8) ≈ 0.612`, so
  reproducing all three is a statement about the branch rather than about one
  point: a solver that merely amplified what it was given could not put the
  maximum in the right place. The residue is the beams' finite temperature and
  moves the right way — widening them at `a = 0.6` gives a monotone approach,
  −7.44%, −5.77%, −4.31%, −3.14%, −2.69%, −2.38% and −2.11% at `vt` from 0.6
  down to 0.15.

  The sharpest assertion is the **stability boundary**, which is qualitative and
  so cannot be laundered by a tolerance: `γ_cold` is exactly zero for `kv₀ ≥ 1`,
  and at `a` = 1.2 and 1.6 the mode decays to 0.052 and 0.000 of its initial
  energy rather than growing. At `a` = 1.0, the cold boundary itself, the warm
  system is still weakly unstable — a factor of 4.9 over `t ≤ 26` — and that is
  asserted as *present* rather than papered over, the boundary being sharp only
  in the cold limit.

  `growth_rate` fits every sample where `damping_rate` fits an envelope, and the
  difference is not an inconsistency: the unstable root here is purely
  imaginary, so the mode grows without oscillating and there are no `log cos²`
  poles — nor, in fact, any local maxima for `damping_rate` to find. Its window
  is set by **amplitude** rather than time, which is what makes it transferable
  across the branch: at `kv₀ = 0.4` a fixed time window gives 10.65%, 6.43%,
  3.99% and 0.19% depending on where it is put, and the amplitude band gives
  2.10% at every `k`. The ceiling also keeps the run inside the solver's
  validity — the field grows with the mode, and a large enough `ε_e` breaks
  `PFC`'s Courant limit in `v`, measured diverging to 1.2e161 before `NaN` at
  `t = 24.1`. Widening the velocity window only postpones that, from `t = 24.1`
  at `±8` to `t = 26.6` at `±16`, which is what identifies the Courant limit
  rather than the boundary as the cause.

  6.7 s, of which 3.8 s is arithmetic.

- **The wakefield study is asserted to run and stay bounded**, and its physics
  now lives in `wakefield` in `test/verification_harness.jl` rather than inline
  in the script. The README has said the example "runs and is stable" for as
  long as it has existed and nothing checked either half — largely because there
  was nothing for a test to call. The script now plots what the shared driver
  returns and reproduces its previous output bit-for-bit.

  The test runs the study at **its own resolution**, in about a second, the
  plotting having been what made the script slow. That is worth insisting on: a
  coarsened proxy is not the same experiment, reporting a peak laser field of
  0.641 against 0.383 and an energy drift of 6.9% against 1.2%, so a test built
  on one would have pinned a different number and called it the study's.

  What it establishes is that the solver runs, stays finite and stays bounded —
  not that the physics is complete, there being no ponderomotive coupling. The
  distinction has teeth in how the assertions divide: the wake and energy
  numbers are bit-identical whether the transverse current is right, wrong by a
  factor of `Δt`, or wrong by thirty-two orders of magnitude, so they constrain
  the longitudinal solver only. **The peak laser field is the single line that
  sees the current**, and is where both documented bugs surfaced — 1.0e22 with
  the `Δt` missing and 44 with the sign flipped, against 0.383 correct. Bounding
  it at 1.0 is the regression net; the rest is a smoke test and says so.

- **The splitting and a real scheme are measured together against an analytic
  answer** (`test/VlasovSolver/test_strang_splitting.jl`). Everything there
  hands the splitting an exact spectral shift, which is the right way to
  isolate the splitting and should stay — but it meant nothing measured what a
  production run actually is. On the same rigid rotation, refining Δx, Δv and Δt
  together over N = 32, 64, 128: cubic `SemiLagrangian` 6.37e-3 → 5.72e-5 at
  orders 3.60, 3.19; `PFC` at 2.26, 2.22; `LaxWendroff` at 1.84, 1.96.

  Each scheme keeps its own spatial order rather than being dragged to the
  splitting's second — the cubic spline reaches 3.19 where Strang alone gives 2.
  That is not the problem being easy: a rigid rotation factors into shears,
  which is what the splitting does, so the commutator error is unusually small
  here and the scheme is what is left.

- **The splitting is reversible, and the round-trip error measures
  dissipation.** Rotation is invariant under `(t, v) → (−t, −v)`, and Strang
  splitting is symmetric, so a forward rotation, a velocity flip, a second
  forward rotation and a second flip must return the initial state. What breaks
  it is the scheme's own dissipation, which has no time-reverse — so this
  separates the schemes far more sharply than the forward error does. At
  N = 64 → 128: cubic `SemiLagrangian` 1.01e-3 → 1.23e-4, `LaxWendroff`
  3.77e-3 → 4.94e-4, `PFC` 1.95e-2 → 3.61e-3, and `Upwind` 4.04e-1 → 2.56e-1.

  Upwind is included **because it fails**: 0.40 is 40% of the peak, it has
  smeared the blob past recovery, and refining barely helps because the
  dissipation is first order. A test every scheme passed would not show that
  this measures dissipation rather than the splitting.

- **The Yee scheme's numerical dispersion relation**
  (`test/MaxwellSolver/test_fdtd_1d.jl`). The update satisfies
  `sin(ωΔt/2) = cfl·sin(kΔx/2)` exactly, and nothing measured it: the suite
  asserted the magic-step case and, at `cfl = 0.8`, only that the deviation from
  a pure translation *exceeded* 0.1 — bounded from below and not from above, so
  a wrong-but-dispersive scheme passed. Measured departure from the closed form
  1.4e-16 to 1.3e-4 over `cfl ∈ {0.5, 0.9, 1.0}` and `kΔx` from 0.016 to 1.41,
  with the physical content quantified: the phase velocity at `λ ≈ 4.4Δx` is
  0.936c at `cfl = 0.5`, 0.981c at 0.9, and exactly c at 1.

  Measured through the mode's own projection `Σ ey·sin(kx)` rather than a point
  sample. A probe at `L/4` reads `sin(mπ/4)`, exactly zero whenever `m` is
  divisible by four — at `m = 20` the "signal" was round-off and the fitted
  frequency came out 4.2× too high.

- **PML reflection against layer thickness.** The existing test measures one
  configuration, which establishes that the layer absorbs but not that the σ
  ramp is what absorbs. Sweeping the thickness at σ_max = 1e3: R = 8.73e-3,
  2.97e-5, 2.16e-8, 3.74e-10, 5.84e-12 at 2, 4, 8, 16 and 32 cells — four orders
  between 2 and 8 cells, which is the cubic ramp working, then a slower fall as
  the limit stops being the layer and becomes the discretisation of the ramp.

- **`BGK` relaxes at the rate it is given.** Both limits of the update were
  pinned bit-for-bit and the moments asserted conserved; the rate in between —
  which is what `τ` means — was not. It holds *exactly*, and the reason ties two
  facts together: `M` is built from the conserved moments, so it is the same
  vector at every step, and `f_k − M = (f_0 − M)·exp(−kΔt/τ)` is then an
  algebraic identity rather than an approximation. A fitted rate that missed
  `1/τ` would mean the moments had moved. Recovered to between 1.4e-10 and
  2.1e-10 at τ = 0.5, 1.0e-13 and 9.4e-13 at τ = 1.0, and 5.1e-15 and 5.8e-14 at
  τ = 2.0 — ranges rather than values, because the fit takes the logarithm of a
  difference that has cancelled to a millionth of its operands, so the residue
  is round-off and tracks the summation order the compiler happens to pick.
  Toggling `--check-bounds=yes` alone moves the middle column by a factor of
  eleven. The assertion sits at `rtol = 1e-8`, fifty times the worst of them.

- **The schemes are compared with each other, not only measured separately**
  (`test/test_comparison.jl`). `benchmark/` times each one in isolation and
  `test_convergence` establishes that each has the order it claims; neither said
  which to reach for, and the answer is not a single name. L² error after one
  traversal at N = 512:

  | scheme | sine | gaussian | square |
  |---|---|---|---|
  | `SemiLagrangian` cubic | 2.46e-8 | 3.55e-6 | 2.01e-2 |
  | `PFC` | 2.30e-7 | 3.30e-5 | 2.65e-2 |
  | `SemiLagrangian` quadratic | 5.02e-6 | 1.37e-4 | 2.72e-2 |
  | `LaxWendroff` | 4.68e-5 | 1.28e-3 | 4.53e-2 |
  | `Godunov`+`VanLeer` | 8.55e-5 | 9.50e-4 | 2.77e-2 |
  | `Upwind` | 8.08e-3 | 4.11e-2 | 6.32e-2 |

  The ranking is strictly ordered by scheme order on the sine, over five
  decades, and **collapses on the discontinuity** to a factor of 3.1 between
  best and worst. The two second-order schemes change places there.
  `LaxWendroff` is linear and so, by Godunov's theorem, cannot be monotone; it
  rings, and falls behind `Godunov`+`VanLeer`, whose limiter buys monotonicity
  with the extrema it clips on the sine. That is the most useful single thing
  to know when choosing a scheme: no ordering of these survives a change of
  problem class. (`Godunov`+`VanLeer` read 7.03e-3, 1.61e-2 and 8.45e-3 while
  its flux lacked a `(1 − |c|)` factor. That put it second from the bottom on
  the sine and first on the pulse, a swing of six decades against the cubic
  spline. Both came from the missing factor, which steepened every slope,
  smooth or not; see Fixed.)

  Errors are asserted; wall-clock is not, for the reasons `runbenchmarks.jl`
  already sets out. The equivalence of `Upwind`, `Godunov(PiecewiseConstant)`
  and linear `SemiLagrangian` shows up here too, agreeing to 2.9e-14 after 1280
  steps — the consequence of the shared amplification factor
  `test_amplification.jl` proves.

- **The collision operators' complexity class is gated.** The one timing
  assertion in the suite, and it is gateable because it measures a class rather
  than a duration: `BGK` doubles when N doubles and `Landau1P` quadruples, so
  the fitted exponents are 1.00 and 1.99 with nothing between them for noise to
  land on. Five trials under `--check-bounds=yes`, which is what `Pkg.test()`
  passes and so the only mode the assertion runs in, gave BGK ratios of 2.00 to
  2.09 and `Landau1P` 3.91 to 4.07, at times from 13.6 us to 6.4 ms. Under the
  `Coverage` job's `--code-coverage=user` every one of those times grows about
  fourteenfold and the exponents do not move — 0.98 to 1.01 and 2.00 to 2.01 —
  because a uniform slowdown cancels in a ratio. That is what makes this
  gateable where a duration is not, and why the 10 us floor `runbenchmarks.jl`
  sets does not apply: nothing here is compared against a number stored on
  another day. It catches an accidental O(N²) in `BGK` — a moment recomputed
  inside the velocity loop, say — which leaves the allocation gate happy and
  every physics assertion passing.

  Samples are drawn against a time budget per size rather than a fixed count, so
  the short measurements the fit is most sensitive to get thousands of them and
  the long ones get five. The testset costs 0.42 s, or 0.7 s under coverage.

  The advection kernels are **not** gated this way, which is a measurement
  rather than an omission: their per-call times are sub-microsecond, and a
  fourfold size increase on an O(N) kernel measured ratios from 3.29 to 8.00
  against the 4 it should give.

- **A work–precision report** (`benchmark/workprecision.jl`, advisory, exits 0).
  Error against the cost of reaching it, per scheme, per problem class, per
  resolution, with the Pareto frontier at the bottom — the schemes no other
  scheme beats on both axes. On smooth data that is `Upwind`, `LaxWendroff`,
  `PFC` and cubic `SemiLagrangian` spanning 0.12 ms to 40 ms and 8.1e-3 to
  2.5e-8, and it is the same four on the pulse. The quadratic spline and
  `Godunov`+`VanLeer` are dominated on every profile, the latter by `PFC`, which
  on the pulse is both more accurate, 2.65e-2 against 2.77e-2, and about half
  the cost. (`Godunov`+`VanLeer` was on the pulse's frontier, and pushed the cubic
  spline off it, while its flux lacked a `(1 − |c|)` factor; see Fixed.
  `Godunov`+`Superbee`, added since, does the same without it: on the pulse it
  reaches 1.58e-2 at 7.2 ms where the cubic spline takes 49 ms for 2.01e-2.)

  Timing goes through `BenchmarkTools.@belapsed` rather than `@elapsed`, which
  is what keeps it from measuring the compiler — the mistake caught in
  `verification/scheme-comparison.jl` a week ago. The report also quantifies its
  own error bar rather than claiming precision it lacks: `Upwind` and
  `Godunov(PiecewiseConstant)` are the same scheme, so any gap between their
  timings is measurement error, and which of them reaches the frontier has
  already changed between two runs of the script.

- **The suite runs against the Julia prerelease, weekly.** A release candidate
  lands weeks before the release, and that gap is the only window in which an
  upstream change that breaks the package — or that the package turns out to be
  relying on by accident — can still be reported upstream and fixed there rather
  than worked around here afterwards. A `Julia prerelease` job now runs the
  default suite at `setup-julia`'s `version: 'pre'`, which resolves to the latest
  RC, beta or alpha, and to the latest stable when no prerelease exists. Between
  cycles the job is therefore a duplicate of the ubuntu `1` matrix entry; that is
  the price of not re-editing a pin every time a cycle opens.

  `continue-on-error`, because the subject under test is Julia and not Vasilek:
  an RC is allowed to be broken, and a red X on the required checks for someone
  else's unreleased bug would train everyone to ignore the checks. The cost is
  silence — a green run sends no failure mail, so a broken RC is visible only in
  the Actions tab.

  The job also brought a `schedule:` trigger to the workflow, at 05:23 UTC on
  Mondays, without which it could not do what it is for: this repository goes a
  fortnight without a push or a PR often enough that an RC could ship, break the
  package and reach its release with the job never having run. The other five
  jobs run on the same cron rather than being gated off it, which is not only
  simpler — the package commits no Manifest, so every run resolves afresh, and
  the weekly run is the only thing that would catch a new `Interpolations` or
  `FFTW` breaking us between one PR and the next.

- **The verification driver takes its schemes as arguments**
  (`test/verification_harness.jl`). `vlasov_poisson` hard-coded
  `PFCNonUniform` on both directions, so the only way to ask what the physics
  costs under a different scheme was to copy the Strang loop — and a copy drifts
  from the one the tests assert against. It now takes `scheme_x` and `scheme_v`,
  defaulting to what every previous caller got, and returns a NamedTuple rather
  than a pair.

  `line_advector` absorbs the API's one asymmetry: `PFCNonUniform` takes a
  displacement while every other scheme takes a Courant number, which only
  exists on a uniform grid. It divides by the spacing for those and **refuses** a
  non-uniform grid rather than picking one of its spacings and being wrong by the
  ratio between them.

- **Mass, momentum, L² and entropy are measured and asserted.** Total energy was
  the only invariant that ever was. Over the k = 0.5 Landau case, 875 steps:
  mass drifts 2.8e-16 and momentum stays at 5.3e-16 on a mass of 25.1 — both
  round-off, both exact conservation laws the discrete scheme also satisfies.
  L² falls 1.0e-5 and entropy rises 7.5e-6, monotonically at every step. Those
  two are *not* conserved and are not asserted as if they were: an exact Vlasov
  flow preserves both, and the drift is numerical dissipation. What is asserted
  is the direction, since a dissipative scheme can only lose L² and gain entropy.

  **The invariants use the cell-width sum `Σ f ΔvΔx`, not `integrate`.** That is
  the quadrature a flux form conserves; the trapezoid halves the two endpoint
  weights, which no conservation law protects. Measured on the same run, the
  trapezoid reports 2.4e-4 of mass drift and 7.3e-4 of momentum against 2.8e-16
  and 5.3e-16 — twelve orders of magnitude of apparent non-conservation that
  belongs entirely to the quadrature. The energy histories keep `integrate`,
  being compared at half-a-percent tolerances where it cannot matter.

- **Landau damping converges under refinement.** Agreement at one resolution
  inside a 3% band can be two errors of opposite sign meeting in the middle.
  Halving Δx, Δv and Δt together — so the Courant number stays at 0.81 and only
  the discretisation moves — gives γ errors of 6.41%, 1.31%, 0.61% and 0.47%.
  The L² dissipation over the same ladder falls by 11.0x, 7.9x and 7.8x against
  the 8x a third-order scheme predicts, which is what identifies the residual
  error in γ: the fitted rate is the physical damping plus the scheme's own,
  which is why every measurement sits on the high side of the analytic value
  rather than scattering about it.

- **A cross-scheme study on a physical observable**
  (`verification/scheme-comparison.jl`, advisory, exits 0). The benchmark suite
  times a bare kernel and `test_convergence` measures an order on a shifted
  sine; neither says what a scheme costs *in the physics*. Ranked by error in
  the Landau damping rate at k = 0.5: `LaxWendroff` 0.64%, cubic
  `SemiLagrangian` 1.07%, `PFC` 1.31%, `Godunov`+`VanLeer` 1.43%, and upwind —
  with `Godunov(PiecewiseConstant)` and linear `SemiLagrangian`, identical to it
  as `test_amplification` requires — at **48.8%**, its own dissipation being two
  orders of magnitude larger than the physical damping it is trying to measure.
  (`Godunov`+`VanLeer` was 3.99% off in the rate, and 16.0% in the frequency
  where it is now 0.54%, while its flux lacked a `(1 − |c|)` factor; see Fixed.)

  The interesting half is the second table. A 1% perturbation cannot rank the
  schemes on positivity at all: every one returns the same `min f = 1.3e-4`,
  which is only the Maxwellian's tail at `v = ±4`. At 50% amplitude the two
  schemes that *lead* the accuracy table are exactly the two that drive `f`
  negative — `LaxWendroff` to −0.094 and cubic `SemiLagrangian` to −0.098,
  against a peak of 0.6. That is Godunov's theorem arriving in the physics, and
  it is why the harness defaults to `PFC` despite it not leading the first table.
  (`Godunov`+`Superbee`, added since, now heads the first table at 0.39% and
  stays positive. It gets there from below, by anti-diffusion rather than
  accuracy; see Added.)

### Fixed

- **The Landau damping rate was fitted with an estimator that its own window
  chose the answer for.** `ε_e ∝ exp(−2γt)·cos²(ωt + φ)`, and the fit ran a
  least squares over `log ε_e` at *every* sample — so it was fitting
  `log cos²`, which has a pole at every null of the oscillation. The window it
  used, `t ∈ [5.9, 29.9]`, began exactly on a minimum. That is the entire
  reason it reported γ = 0.1498, 2.3% below the tabulated 0.15336; moving the
  start one step, to 6.0, gives 0.1532 on the same data.

  Sweeping plausible windows moved the old estimate over 0.14837 to 0.15755, a
  spread of 6.2% — larger than the 5% tolerance it was being held to, so the
  test was passing on the strength of where its window happened to land.
  Fitting through the local maxima instead removes the `cos²` entirely: the
  same sweep now gives 0.15451 to 0.15571, a spread of 0.8%, consistently about
  1% above the analytic value. That residue is numerical damping and does not
  move when Δx, Δv and Δt are all halved.

  A window-sensitivity assertion is now part of the test, so this class of
  failure cannot come back silently.

- **`local_extrema` skipped the last interior point it was asked about.** The
  helper iterated `2:length(y)-2`, but the strict-interior comparison is
  well-defined up to `length(y)-1`, where the right neighbour is `y[end]`. On
  its own documented terms — "the strict interior local maxima … or minima of
  `y`" — it was one short: `[10, 20, 5]` returned nothing at all, and
  `[1, 3, 2, 4, 1]` returned only the first of its two peaks. No measured value
  moves, because both callers read `ε_e` from `vlasov_poisson`, which pads
  `ε_e[end] = ε_e[end-1]`, so a strict inequality at index `end-1` compares
  against an equal neighbour and cannot fire — an invariant that lived in the
  caller and was not stated on the helper.

### Added

- **Landau damping at three wavenumbers, and the real frequency**
  (`test/test_verification.jl`). One `k` with one fitted `γ` is a single point
  on a curve and cannot separate a solver that reproduces the dispersion
  relation from one that lands near a value at one wavenumber. Both roots are
  now measured at `k = 0.3, 0.4, 0.5`:

  | k | γ | vs tabulated | ω_r | vs tabulated |
  |---|---|---|---|---|
  | 0.3 | 0.01257 | 0.42% | 1.15696 | 0.25% |
  | 0.4 | 0.06646 | 0.51% | 1.28042 | 0.36% |
  | 0.5 | 0.15558 | 1.45% | 1.41372 | 0.14% |

  The real frequency is about four times the sharper of the two — it comes from
  counting nulls, where `γ` comes from fitting an amplitude that numerical
  dissipation also acts on — so it is held to 1% where `γ` is held to 3%.

  Two constraints are documented in the test because they are not obvious and
  cost real time to rediscover. The velocity window must contain the resonance
  at `v = ω_r/k` (2.83, 3.21, 3.87), since that is where the damping comes
  from; and the window then fixes `Δt`, because the fastest row runs at
  `max|v|·Δt/Δx` against `PFCNonUniform`'s Courant limit of 1 — so reaching a
  resonance costs a smaller time step, not just more velocity points. The
  existing `k = 0.5` case ran at 1.019, marginally over; it now runs at 0.815.
  Separately, the fitting window has to stop before the mode reaches the floor
  where recurrence and round-off take over, and that arrives earlier *in units
  of e-foldings* the weaker the damping is: at `k = 0.3` the mode has decayed
  only threefold by `t = 50`, and fitting past it reports γ = 0.0099, 22% low.

- **The frequency estimator is now held to its window as well.** `γ` gained a
  window-sensitivity assertion above; `ω` rested on a counting assumption of its
  own and had none. Its spacing is averaged between the first and last minimum
  over `length(m) − 1` intervals, so one null missed on a near-tie — or one
  spurious null off numerical noise — rescales the answer with nothing to say
  so, and agreement with the analytic value on a single window could be luck.
  Both windows are now read for `ω` too: the three Landau cases agree to 0.05%,
  0.15% and 0.10% against a 1% threshold, and the Bohm–Gross run to 0.006%
  against 0.2%. It costs no extra solve, only a second reading of the same
  `ε_e`.

  The thresholds sit in a gap worth naming. Below them is the estimator's own
  floor: a null is located only to within `Δt`, so two windows disagree by about
  `Δt/span` whatever the physics does — 0.4% at `k = 0.5`, and 0.05% for the
  long Bohm–Gross window. Above them is the failure being tested for: losing one
  null of ten rescales the spacing by 11%, one of sixty by 1.7%.

- **The plasma oscillation frequency is asserted against Bohm–Gross.**
  `docs/normalization.md` states that the plasma-oscillation study verifies the
  analytic plasma frequency; nothing in the repository measured a frequency
  anywhere, only the energy drift — which a solver oscillating at entirely the
  wrong rate passes without difficulty. Measured ω = 1.005719 against
  `√(1 + 3k²) = 1.005904`, 0.018%, where the cold `ωₚ = 1` is 0.57% away. The
  factor of thirty is the point: the test fails if the thermal correction is
  dropped rather than merely preferring it to be there. Costs one second at
  `t ≤ 200`, where the energy test needs 3000 and forty times as long.

- **The extended verification suite now runs in CI.** `VASILEK_EXTENDED=1` was
  named in the README, the changelog and three test files, and set by none of
  the four CI jobs — so the body of `test/test_verification.jl` had never
  executed on a runner, and the three claims in the README's verification table
  held only when someone remembered to export the variable by hand. There is now
  an `Extended verification` job that sets it. Per-PR rather than nightly: the
  benchmark suite is advisory because wall-clock timing on a shared runner is
  unreliable at any tolerance worth having, and that argument does not carry over
  to deterministic numerics. Measured 2m08s for the whole job against 51s for the
  default suite; the three assertions pass.

- **Von Neumann amplification factors, against closed forms**
  (`test/test_amplification.jl`). Each scheme is linear, so a Fourier mode is an
  eigenvector and one step multiplies it by a `g(kΔx, c)` available in closed
  form. Six schemes are asserted mode by mode against symbols derived from the
  update formulas: `Upwind`, `Godunov(PiecewiseConstant)` and
  `SemiLagrangian(LinearSpline)` share `1 − c(1 − e^{−iθ})`; `LaxWendroff` and
  `Godunov(PiecewiseLinear, NoLimiter)` share `1 − ic·sinθ + c²(cosθ − 1)`; and
  `PFC`'s unlimited branch is its third-order flux. Measured agreement over
  `m ∈ {1,2,4,8,16,24,31}` and `c ∈ {0.4, 0.8}`: 1.6e-14 to 1.8e-14 for the
  three-point schemes, 4.6e-14 for `PFC`. (The first was quoted as 8.9e-15;
  Julia 1.10.12 and 1.13.0 both give 1.8e-14 for upwind.)

  This is sharper than what guarded these kernels before. `test_convergence`
  fits a slope to ±0.15 and `test_golden` pins one dataset for eight steps;
  this pins every mode from the fundamental to the grid scale against an
  analytic value, and because the comparison runs elementwise it asserts at the
  same time that the output *is* a pure mode. One consequence worth naming: the
  three claimed scheme equivalences now hold mode by mode rather than on one
  profile. (There was a fourth symbol, the centred `1 − ic·sinθ` of
  `Godunov(PiecewiseLinear, NoLimiter)`, with `|g| > 1` for every mode — worst
  1.077 per step at `c = 0.4`. That stated its unconditional instability
  analytically, until its flux gained the `(1 − |c|)` factor that makes it
  `LaxWendroff`'s; see Fixed.)

- **Dissipation and dispersion, per scheme per wavenumber**, from the same
  symbols at no extra cost: `|g|` is the amplitude lost per step and
  `arg(g)/(−cθ)` the relative phase velocity. The two order the schemes
  *differently*, which is the point. At `c = 0.4`, `λ = 64Δx`: upwind loses
  1.16e-3 per step where `LaxWendroff` and `PFC` lose 2e-6. At `λ = 4Δx` the
  phase error is 2.4% for `PFC` against 29% for `LaxWendroff` — `LaxWendroff`
  keeps a marginal mode's amplitude and puts it in the wrong place, where upwind
  removes it instead. Also asserted from the symbol alone: the amplitude retained
  after one full domain traversal, 0.831 for upwind against 0.9998 for
  `LaxWendroff` and `PFC`, which is the same 17% loss that shows up in the
  advection error tables, arrived at independently.

- **Free streaming: phase mixing and recurrence**
  (`test/VlasovSolver/test_free_streaming.jl`). With no field the Vlasov
  equation has an exact solution, and a spatially modulated Maxwellian gives a
  density mode decaying as `exp(−k²t²/2)` — a Gaussian, not an exponential,
  because the mode is sheared into fine velocity structure rather than
  dissipated. On a discrete velocity grid the sheared modes rephase and the mode
  returns at `T_R = 2π/(kΔv)`, a number set by `Δv` alone that no scheme can
  move.

  This is the only exact kinetic solution available without a field solve, and
  it exercises the x-sweep together with the velocity moment outside the
  extended gate — that path previously ran only inside `verification_harness.jl`.
  Measured relative error in the mode amplitude over `t ≤ 4`: 2.7e-6 for cubic
  semi-Lagrangian, 3.3e-5 for `PFC`, 1.2e-3 for `LaxWendroff`, 1.7e-2 for
  upwind. Four orders of magnitude on a quantity with physical meaning, and the
  ordering is asserted rather than only printed, errors being deterministic
  where timings are not. Recurrence is measured on a deliberately coarse
  velocity grid, which brings `T_R` from 125.7 to 31.4 and makes the test four
  times cheaper without weakening it: both schemes peak at `t = 31.42` against
  `T_R = 31.4159`, recovering 1.0000 and 0.9991 of the initial amplitude.

  The Courant limit binds here in a way it does not in the one-dimensional
  advection tests: the fastest velocity row carries `max|v|·Δt/Δx`, so the
  velocity window and the time step are not independent. Everything runs at
  0.611.

- **The Yee leapfrog's conserved energy is asserted.** `E` sits at integer steps
  and `H` at half-integer ones, so the conserved quadratic form is staggered in
  time too — `‖E^{n+1}‖² + ⟨H^{n+1/2}, H^{n+3/2}⟩`, equivalently
  `⟨E^n, E^{n+1}⟩ + ‖H^{n+1/2}‖²`. Measured relative drift over 2000 steps:
  3e-15 at cfl = 0.5, 0.8 and 0.99. The naive `‖E‖² + ‖H‖²` at a single instant
  is **not** conserved — it ranges over 12% to 24% across the same runs — and
  the test asserts that too, so nobody reaches for it later.
- **PML absorption is measured rather than eyeballed.** The previous test
  checked the field ended up near zero, which a badly-but-symmetrically
  absorbing layer also passes. Now: reflection coefficient 6.25e-9 rightgoing
  and 6.26e-9 leftgoing, asymmetry 1.0010 — the two ends have different index
  arithmetic, so an off-by-one in one of them was invisible to a test that only
  looked one way — and a no-PML control leaving 0.999 of the pulse on the grid,
  so the numbers are the layer working rather than the pulse having left.
- **`BGK` satisfies the H-theorem**, to machine precision once the velocity
  window resolves the relaxed state. The most negative single-step increment
  over 200 steps is -1.2e-4 at ±6, -3.1e-12 at ±10 and -4.4e-16 at ±14, while
  refining Δv does not move it — which is what identifies the truncation rather
  than the operator as the cause, and the test asserts both halves of that.
- **`∂f∂v` is tested directly**: second order in the interior, first at the
  ends, and exact on a linear profile including on a non-uniform grid. It was
  factored out of `Landau1P` after the second copy was found differentiating at
  the wrong index, and until now was only ever exercised through the operator
  that misused it.
- `YeeMesh1D` shape and element type, and a `Bounds checking` CI job running the
  suite under `--check-bounds=yes` so the `@inbounds` on the advection kernels
  is non-binding for one run. It passes today; the point is that nothing else
  could have told us.
- **The README's usage example is executed by the suite** and so cannot drift.
  It had: the block referred to `src`, `dest` and `courant` without defining any
  of them, and the file carried two near-identical `## Usage` sections.

- **Reentrancy is asserted** (`test/test_threading.jl`), and a `Multithreaded`
  CI job runs the suite at `JULIA_NUM_THREADS=4` so it means something. Sweeping
  many independent lines with one shared scheme value and a per-task workspace
  must reproduce the serial result *bit-for-bit*, in both Courant directions,
  for every scheme and for `BGK`. This is the property the 0.2 refactor existed
  for — the old closures captured a single shared scratch buffer, which is why
  the parallel 2D2P goal was unreachable — and it was the one claim with no test.
  Workspaces are allocated inside the task rather than indexed by
  `Threads.threadid()`: tasks migrate, and `maxthreadid()` exceeds `nthreads()`
  whenever an interactive pool exists (8 against 4 here), so a `threadid()`-keyed
  pool is both a race and an out-of-bounds read waiting to happen. The test is
  the pattern the 2D2P sweeps should copy.
- **`BGK` moment conservation** (`test/BoltzmannSolver/test_damping_1v.jl`).
  Leaving `n`, `u` and `T` alone is the entire content of the operator and
  nothing checked it; worse, every existing test started from data symmetric in
  `v`, so `u` was zero throughout and the mean-velocity computation was never
  exercised at all. Measured over 100 steps on `v ∈ [-10, 10]`: machine
  precision for skewed and drifting initial data, 2.5e-10 for a bi-Maxwellian.
  Adds an arbitrary Maxwellian as a fixed point at four `(u, T)` pairs, and both
  limits of the update — `Δt ≪ τ`, and `Δt ≫ τ` giving the local Maxwellian
  bit-for-bit — plus the intermediate `f·e + (1−e)·M`, also bit-for-bit.
  The residual at `T = 2` is quantified as the velocity window rather than left
  as a footnote: ±8 gives 1.6e-5 where ±10 gives 3.7e-9, which is the tail the
  trapezoid loses, and it sets the window callers need.
- **Contract tests** for `advect!` (`test/test_contracts.jl`): the argument
  checks below, workspaces carrying no state between calls, `@inferred` on every
  kernel, and element-type behaviour — including the pinned fact that
  `workspace(::SemiLagrangian, n)` hard-codes `Vector{Float64}`, so a `Float32`
  line is prefiltered in double precision.

- **Direction-symmetry suite** (`test/test_symmetry.jl`). Every other advection
  suite ran at `c = 0.4 > 0`, and `c > 0` and `c < 0` are separate branches in
  every scheme that has one; only a single-step check ever touched the negative
  one. Instead of duplicating four suites, this asserts the identity that ties
  the branches together, `R∘A(+c) == A(−c)∘R` — exact for `Upwind` and
  `Godunov`, 4.4e-16 for `LaxWendroff` and `PFC`, 4.5e-14 for the
  semi-Lagrangian family, measured over five datasets × four Courant numbers.
  That it suffices was checked rather than assumed: backward-direction
  convergence orders and invariants agree with the forward ones to three
  decimals. Also covers constant-preservation, `c = 0` as the identity (which
  runs the negative branch, `c > 0` being false at zero), periodicity of
  `SemiLagrangian` in the Courant number — the only test of its periodic
  extrapolation path — and its stability at `c = 3.7`, where every other scheme
  diverges.
- **`StrangSplitting` has tests** (`test/VlasovSolver/test_strang_splitting.jl`).
  It was executed zero times in a default run, being reachable only through
  `verification_harness.jl` behind `VASILEK_EXTENDED=1`. Second order in Δt is
  asserted on rigid rotation, measured 2.003/2.001/2.002 against both a fine
  reference and the analytic answer. The advection operator there is an exact
  spectral shift, deliberately: a scheme's spatial error does not vanish as
  Δt → 0 — semi-Lagrangian interpolation error accumulates with the step count —
  so measuring the order against a real scheme measures the scheme. Also pins
  the non-square `Nx ≠ Nv` transpose bookkeeping, and the postcondition that
  **`f[2]` lags `f[1]` by the final half-step on return**, which every caller in
  this repository reads diagnostics off.
- **`PFC` bounds check is tested** (`test/test_contracts.jl`). The `checked`
  type parameter had a paragraph of docstring, the memory of a 740× bug, and no
  test. Asserts that the check fires in both branches and at neither end of the
  closed interval, that `checked = false` is bit-for-bit identical on valid data,
  and that the element type survives `promote` — including that integer bounds
  *widen* a `Float32` pairing rather than narrowing to it.
- **The Poisson solver's exact discrete identity is asserted**, closing a gap
  between `docs/normalization.md` — which already claimed the test asserted it
  "to machine precision" — and the test, which compared one mode at
  `atol = 1e-3`. `E_num == sinc(kΔx)·E_exact` now holds to 1.2e-13 across even
  and odd `N`, six modes and five grid spacings, including one mode below
  Nyquist. The departure is normalised by the continuum amplitude `1/k` rather
  than by the output: near Nyquist the sinc factor suppresses the field itself,
  and dividing by the output would report 7.3e-12 for the same 2e-15 absolute
  error the solver has everywhere. Adds linearity, invariance to a net charge,
  bit-for-bit repeatability across calls, and second order in Δx.
- **FDTD at the magic time step.** At `cfl = 1` the 1D Yee update has no
  numerical dispersion and translates a pulse exactly: measured 1.1e-19 after
  20 steps against a plain `circshift`, where a mis-signed or mis-indexed term
  would leave an O(1) residue. `cfl = 0.8` is run beside it, off by 0.83, so the
  test visibly has teeth. The Courant limit is asserted too — bounded at
  `cfl = 0.99` and 1, non-finite at 1.02 and 1.2 — which nothing covered.

- **Extended verification** (`test/test_verification.jl`), gated behind
  `VASILEK_EXTENDED=1`. The claims the verification documents made in prose
  are now assertions: Landau damping at k = 0.5 within 5% of the tabulated
  0.1533 (measured 0.1498), and energy drift below 0.5% on the uniform grid
  and 6% on the non-uniform one at t = 3000 (measured 0.38% and 4.85%).
- [docs/normalization.md](docs/normalization.md): the unit conventions, the
  wavenumber convention that the Poisson bug came from, the FDTD current
  convention that the wakefield instability came from, and why the PFC bounds
  have no default.
- `PFC(; fmin, fmax, checked = false)` compiles the bounds check away. The
  check moved into `advect!` when the scheme stopped seeing data at
  construction, costing a minimum/maximum pass per call -- 13% of the step at
  N = 10000, which I had not measured at the time. Now a type parameter, on
  by default.

### Changed

- **`FDTD1D` drives only the interior nodes, and the PEC boundaries are
  written down.** Both ends of the grid are perfect electric conductors: no
  update touches `ey[1]`, `ez[1]`, `ey[end]` or `ez[end]`, so they hold the zero
  `YeeMesh1D` gives them. That was already true of the curl loops, but only by
  omission and nowhere stated — `update_ey!` and `update_ez!` added the current
  over `1:Nx`, so node 1 was driven while node `Nx+1` was not.

  Injecting a current into a perfect conductor is meaningless, and doing it at
  one end only broke the symmetry the conserved staggered energy rests on: the
  discrete curls are adjoint precisely because the boundary terms vanish. The
  current is now applied over `2:Nx`, exactly the set of dynamic nodes. The
  arrays still span all `N+1` nodes so they index alongside `f.ey`; their first
  and last entries are ignored.

  Callers should seed the interior: writing into an end node does not launch a
  wave, it changes the boundary condition to `E = const`. The convention is now
  in [docs/normalization.md](docs/normalization.md) and in the `YeeMesh1D` and
  `make_advance_fields` docstrings, and the tests assert the observable
  consequence rather than the loop bounds — a pulse reaching either wall returns
  inverted, with a measured reflection coefficient of −0.9998.

  No effect on the verification runs: `wakefield.jl` is the only caller passing
  a nonzero current, and its density profile vanishes at both walls.

- **`advect!` validates its arguments.** Three ways of calling it wrongly used
  to produce a plausible wrong answer rather than an error, which is the class
  of failure this package has already paid for twice — the two `PFC` overloads
  that disagreed by a factor of 740, and the `LaxWendroff` methods that sized
  their loops from a captured array.

  * `dest === src` is rejected. Every scheme except `SemiLagrangian` and
    `PFCNonUniform` reads neighbours of `src` that an aliased `dest` has already
    overwritten — measured 0.013 off for `Upwind`, 0.21 for `PFC`. The two that
    survived did so by accident of holding a scratch copy, not by contract, so
    aliasing is rejected uniformly rather than left scheme-dependent.
  * `length(dest) != length(src)` is rejected. It used to be half-loud:
    `length(dest) < length(src)` truncated silently in the finite-difference
    schemes, wrapping periodically at the shorter length, while the other
    direction threw `BoundsError`.
  * A workspace of the wrong size is rejected. Undersized already threw;
    **oversized did not**, and that was the dangerous one — `interpolate!`
    prefilters the whole buffer, so a `SemiLagrangian` handed a workspace built
    for a longer line returned garbage (max|Δ| ≈ 0.4) without complaint. The
    check is exact equality, so a single large workspace can no longer be shared
    across lines of differing length.
  * Fewer than three cells is rejected. `Godunov(PiecewiseLinear())` reaches two
    neighbours either side and threw `BoundsError` at `n = 2` while every other
    scheme quietly returned a degenerate answer; the floor is now uniform.

  The checks are three comparisons outside the loop with `@noinline` error
  paths. The allocation gate still reports zero bytes for every kernel.

- **The verification notebooks are runnable scripts.** `.jmd` + Weave becomes
  Literate-format `.jl` that executes directly and writes its own figures.
  Removes the Weave dependency and 330 kB of HTML generated in 2021 against
  code that is no longer in the package.
- **`wakefield.jl` produces output and is stable.** It used to plot the final
  snapshot where a time-space matrix was wanted, which never raised because a
  script never renders. Making it visible showed a peak field of 1.0e22
  against a source amplitude of 1.0; a missing Δt on the current and a sign
  taken from the verified longitudinal convention bring it to 0.38.

  It still has **no ponderomotive coupling** -- the laser does not drive the
  wake, evidenced by the wake coming out bit-identical whether the transverse
  current is right or wrong by 32 orders of magnitude. Left to the author.

### Breaking

- **Advection schemes are types, and `generate_solver` is gone.** `advect!`
  dispatches on an immutable scheme value; scratch memory, where a scheme
  needs any, comes from `workspace` and is passed explicitly. See
  [docs/migration-0.2.md](docs/migration-0.2.md).

  No compatibility shim, contrary to the original plan: the module names
  `generate_solver` lived under are the new type names, so the two APIs
  cannot coexist in one namespace.

  The reason is reentrancy. A 2D2P step sweeps O(N²) independent lines per
  direction, which is work for `Threads.@threads`; the old closures captured
  a single shared scratch buffer, so they could not be threaded over them.
  The parallel 2D2P goal in the README was unreachable without this.

- The `Symbol`-valued options are singleton types: `:Riemann_linear` becomes
  `PiecewiseLinear()`, `:VanLeer` becomes `VanLeer()`, `:Cubic` becomes
  `CubicSpline()`. A typo is now a `MethodError` where it is written.
- **Collision operators follow the same shape**: `BGK(τ)` and
  `collide!(dest, src, op, v, Δt, ws)` replace `BGK.generate_solver`, which
  mutated its argument in place. `Landau1P` is unexported and experimental.
- `Limiters` is gone; the limiters are callable types in `Advection`.

  Numerics are unchanged. Every scheme is bit-for-bit identical to 0.1 in
  both directions, and the golden values recorded before the change still
  match. (`Godunov(PiecewiseLinear())` has since changed on purpose, gaining the
  `(1 − |c|)` factor its flux lacked in 0.1; see Fixed.)


### Added

- **Allocation gate** (`test/test_allocations.jl`): blocking, deterministic, and
  asserting zero for every kernel that achieves it. `@allocated` does not flap
  the way wall-clock timing does, so it can gate CI.
- `--strict` and `--rebaseline` flags for the benchmark driver, and
  `PFCNonUniform` and `FDTD1D` added to the benchmark suite.

- **Golden-value regression tests** (`test/test_golden.jl`): bit-for-bit
  reference output for all nine advection scheme/option combinations, eight
  steps each. Convergence rates are far too loose to catch an index swapped
  by one; bit-identity is not. Regenerate with `test/generate_golden.jl`, and
  only when a numerical change is intended. Advection-only by design — the
  collision operators go through libm and already differ between Julia
  1.10 and 1.12.
- **Order-of-accuracy suite** (`test/test_convergence.jl`): global convergence
  order for eight schemes, measured before being asserted. Replaces reliance
  on magic `atol` values at a single resolution.
- **Structural invariants** (`test/test_invariants.jl`): mass conservation,
  total variation, PFC positivity and the maximum principle, and the Courant
  limit. None of this was covered.

- Continuous integration on GitHub Actions: tests on Julia `lts` and `1` across
  Linux, Windows and macOS, plus a coverage job reporting to Codecov.
  Verified locally on both ends of the matrix (1.10.12 and 1.12.7).
- Dependabot updates for GitHub Actions.
- `Aqua.test_all` in the test suite.
- Separate `benchmark/Project.toml` and `verification/Project.toml`.
- Advection tests now assert that both `generate_solver` overloads produce
  bit-identical output — the invariant `PFC` had been violating.

### Fixed

- **The `(h, h₀)` methods of `LaxWendroff` sized their loops from the captured
  array** rather than from the pair passed to them, so a differently sized
  pair read and wrote out of range. The same mistake `PFC` had. Both now
  have a regression test, which neither did.
- `Godunov.generate_Φ` and `Limiters.generate_limiter` returned `nothing` for
  an unrecognised option, surfacing as a call on `nothing` from inside the
  hot loop. Both now throw `ArgumentError` at construction.

- **`SemiLagrangian` rebuilt its interpolation on every step.** `interpolate`
  copies its input and runs the B-spline prefilter on the copy, costing 80 kB
  per step at N = 1000 for linear and 173 kB for quadratic and cubic. The
  buffer is now allocated once and prefiltered in place: linear drops to zero,
  quadratic and cubic to about 100 kB, all of it inside the periodic prefilter
  in Interpolations itself. Output is bit-identical, confirmed by the golden
  tests.
- **`BGK` allocated six temporaries per step** (39 kB at 801 velocity points).
  Two reusable buffers and fused in-place broadcasts take it to zero, with
  bit-identical output over 100 steps.
- Benchmark measurements were invalid: `@benchmarkable` interpolated only `N`,
  leaving a non-`const` global `Dict` lookup inside the timed region. The
  closures are now interpolated directly.
- The benchmark driver ran only the Vlasov group, so the Maxwell group was
  tuned and never executed. `results.json` stored every raw sample (4.4 MB);
  it now stores the estimate (13 kB).
- The Cyrillic `c` (U+0441) in the benchmark keys is gone.

- **`FDTD1D` read `f.hz[i-1]` in an `hy`-only stencil for `ez`.** Because `hz`
  is non-zero throughout the existing tests, this actively injected `ez`,
  which then drove `hy`, which fed back: seeding only the y polarization at
  amplitude 1.0 gave `max|ez| = 7.84` after ten steps. The `(ey, hz)` pair
  survived only because it never reads `hy`. Numerically breaking for any
  z-polarized run.
- `FDTD1D`'s `PML` default argument passed `Δt` and `Δx` in the wrong order.
  Added a keyword constructor so the mistake is unrepresentable.
- **`Landau1P` differentiated at the wrong index.** The inner loop computed
  `Δfⱼ` but branched on `i` and read `f₀[i±1]`, so `Δfⱼ ≡ Δfᵢ`; and `J[j]` was
  never assigned when `i == j`, leaving a stale or uninitialised value that
  was then integrated. `L` and `Tₜ` are now keyword arguments.
- **`PFCNonUniform` computed its limiter from the global spacing ratio,** so
  one refined region tightened it everywhere (ξ ≈ 0.571 against 2 on the
  notebooks' velocity grid). Now per cell triple. Energy drift of the
  non-uniform plasma-oscillation case falls from a reported 12% to 4.85% at
  t = 3000; the uniform grid is bit-identical, as it must be.
- `PFCNonUniform` hardcoded `fₘᵢₙ = 0`, `fₘₐₓ = 1`; now required keywords, as
  in `PFC`. **Breaking API change.**
- **`PoissonFourier1D` returned a field too large by `1/Δx²`.** `rfftfreq(n)`
  yields cycles per *sample*, so `Δx` never entered the spectrum while the
  φ→E step divided by it once. Numerically breaking: any result computed with
  this module before this change was wrong by that factor. Neither verification
  notebook was affected — both define their own Poisson solver.
- **`PFC`'s two `generate_solver` overloads computed different things.** The
  keyword form defaulted `fₘₐₓ = 1` while the Courant-baked form used
  `maximum(f₀)`; on data exceeding 1 the limiter returned a spurious slope
  (measured divergence: a factor of 740). `fₘᵢₙ` and `fₘₐₓ` are now required
  keywords on both, asserted to bracket `f₀`. **Breaking API change.**
- `PFC`'s three-argument `solve!(h, h₀, c)` sized its loops from the captured
  array rather than from `h`, breaking whenever the two differed.
- The Poisson test sampled 1001 points of a period-10 sine, so the input was
  not periodic on the DFT window and leaked spectrally.
- `test_damping_1v.jl` rebound a local in its stepping loop, so the closure kept
  reading the original buffer and the "100 time steps" repeated step one a
  hundred times. The `Landau1P` case is now `@test_broken`: with the loop
  actually running, the known index bug surfaces (deviation 0.0022 → 0.31).
- `test_fdtd_1d_current` had its assertion commented out and only printed.
  Restored; it passes.

### Changed

- `Upwind` and `LaxWendroff` reduced to a single numerical kernel each, with
  both `generate_solver` methods as wrappers. The kernel had been written out
  twice per module, differing only in whether the Courant number was captured
  or passed. Output is bit-identical, confirmed by the golden tests.

- Tests load the package (`using Vasilek`) instead of including source files by
  relative path and dispatching through `eval`. `src/Vasilek.jl` was previously
  never exercised by its own test suite.
- Runtime dependencies cut from nine to four: `FFTW`, `Interpolations`,
  `LinearAlgebra`, `NumericalIntegration`. `DSP` and `Weave` were imported
  nowhere; `Plots` only inside dead branches; `BenchmarkTools` belongs to
  benchmarks. Added the missing `[compat]` section and `[extras]`/`[targets]`.
- `@fastmath` removed from the advection kernels. Measured at N = 10000 it
  bought nothing (1.55→1.51 µs, 2.47→2.38 µs, 48.2→48.1 µs) and output is
  bit-identical, but it licensed LLVM to reassociate floating-point
  arithmetic — and bit-for-bit reproducibility is what the planned
  golden-value tests depend on. `@inbounds @simd` are kept.
- `Limiters` hoisted to the top level; `Godunov` now uses `using ..Limiters`
  rather than nested-including it, so the test process no longer holds two
  distinct `Limiters` modules.
- `StrangSplitting` added to the package include list. It is one of the three
  features the README advertises and both verification notebooks depend on it,
  yet it was reachable from neither.
- `Manifest.toml` is gitignored rather than tracked, as it should be for a
  library. The previously committed manifest was still in v1.0 format despite
  the "upgrade to Julia v.1.11" commit, so it never recorded which Julia had
  resolved it.

### Removed

- Five files containing a single no-op function each (`Maxwell1D`, `Maxwell2D`,
  `Poisson2D`, `Advection2D`, `Landau2P`) and the three `src/Mesh` types — bare
  structs with untyped fields, referenced by nothing and absent from the include
  list. A typed grid abstraction belongs with the 2D work.
- The `plot_needed` / `plot!` machinery in tests. Never exercised, and the only
  reason the test environment needed Plots and GR.

### Known issues

- `Landau1P`'s closure is inconsistent: the numerator carries the transversal
  estimate `2Tₜ` while the denominator uses the longitudinal `|vᵢ-vⱼ|³`, leaving
  a non-integrable singularity at `i ≈ j`. The collision integral grows under
  grid refinement instead of converging (0.0023, 0.0053, 0.0080, 0.0104 at
  Δv = 0.4, 0.2, 0.1, 0.05) and a Maxwellian is not a stationary point. Two
  `@test_broken` assertions record this. A consistent closure would use
  `(Δv² + 2Tₜ)^(3/2)`; that is a physics decision, not a coding fix.
- `Landau1P` differences a cell-centred `I` rather than staggered fluxes, so
  mass is not conserved to machine precision. Measured drift over 100 steps is
  2e-10, which the test asserts as a bound.
- Coverage is collected but discarded: `CODECOV_TOKEN` is not configured.
