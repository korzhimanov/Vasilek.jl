# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project has not been released; entries below describe work on `master`.

## [0.2.0] - unreleased

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

  Measured at three wavenumbers, with the beams at `vt = 0.3`: γ = 0.30244,
  0.34229 and 0.31232 against 0.30819, 0.35339 and 0.31134 — 1.87%, 3.14% and
  0.31%. **`γ(a)` is non-monotone**, peaking at `a = √(3/8) ≈ 0.612`, so
  reproducing all three is a statement about the branch rather than about one
  point: a solver that merely amplified what it was given could not put the
  maximum in the right place. The residue is the beams' finite temperature and
  moves the right way — widening them to `vt` = 0.6, 0.5, 0.4, 0.3 gives 3.53%,
  2.03%, 0.92% and 0.11% on a fixed window.

  The sharpest assertion is the **stability boundary**, which is qualitative and
  so cannot be laundered by a tolerance: `γ_cold` is exactly zero for `kv₀ ≥ 1`,
  and at `a` = 1.2 and 1.6 the mode decays to 0.053 and 0.000 of its initial
  energy rather than growing. At `a` = 1.0, the cold boundary itself, the warm
  system is still weakly unstable — a factor of 4.9 over `t ≤ 26` — and that is
  asserted as *present* rather than papered over, the boundary being sharp only
  in the cold limit.

  `growth_rate` fits every sample where `damping_rate` fits an envelope, and the
  difference is not an inconsistency: the unstable root here is purely
  imaginary, so the mode grows without oscillating and there are no `log cos²`
  poles — nor, in fact, any local maxima for `damping_rate` to find. Its window
  is set by **amplitude** rather than time, which is what makes it transferable
  across the branch: at `kv₀ = 0.4` a fixed time window gives 9.76%, 5.63%,
  4.35% and 0.56% depending on where it is put, and the amplitude band gives
  1.87% at every `k`. The ceiling also keeps the run inside the solver's
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
  | `Godunov`+`VanLeer` | 7.03e-3 | 1.61e-2 | **8.45e-3** |
  | `Upwind` | 8.08e-3 | 4.11e-2 | 6.32e-2 |

  The ranking is strictly ordered by scheme order on the sine, over five
  decades — and **inverts on the discontinuity**. `Godunov`+`VanLeer` sits
  second from the bottom on smooth data, because the limiter clips smooth
  extrema, and first on the pulse; against the cubic spline it goes from 2.9e5
  times worse to 0.42 times as bad, a swing of six decades in relative standing.
  That is Godunov's theorem as a measurement, and it is the most useful single
  thing to know when choosing a scheme: no ordering of these survives a change
  of problem class.

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
  2.5e-8; on the pulse `Godunov`+`VanLeer` joins it and the cubic spline drops
  off, being both slower and less accurate there than `Godunov`+`VanLeer`. The
  quadratic spline is dominated on every profile.

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
  mass drifts 2.8e-16 and momentum stays at 1.7e-15 on a mass of 25.1 — both
  round-off, both exact conservation laws the discrete scheme also satisfies.
  L² falls 9.8e-6 and entropy rises 7.4e-6, monotonically at every step. Those
  two are *not* conserved and are not asserted as if they were: an exact Vlasov
  flow preserves both, and the drift is numerical dissipation. What is asserted
  is the direction, since a dissipative scheme can only lose L² and gain entropy.

  **The invariants use the cell-width sum `Σ f ΔvΔx`, not `integrate`.** That is
  the quadrature a flux form conserves; the trapezoid halves the two endpoint
  weights, which no conservation law protects. Measured on the same run, the
  trapezoid reports 1.6e-4 of mass drift and 7.3e-4 of momentum against 2.8e-16
  and 1.7e-15 — three orders of magnitude of apparent non-conservation that
  belongs entirely to the quadrature. The energy histories keep `integrate`,
  being compared at half-a-percent tolerances where it cannot matter.

- **Landau damping converges under refinement.** Agreement at one resolution
  inside a 3% band can be two errors of opposite sign meeting in the middle.
  Halving Δx, Δv and Δt together — so the Courant number stays at 0.81 and only
  the discretisation moves — gives γ errors of 5.83%, 1.45%, 0.62% and 0.47%.
  The L² dissipation over the same ladder falls by 10.2x, 7.8x and 7.8x against
  the 8x a third-order scheme predicts, which is what identifies the residual
  error in γ: the fitted rate is the physical damping plus the scheme's own,
  which is why every measurement sits on the high side of the analytic value
  rather than scattering about it.

- **A cross-scheme study on a physical observable**
  (`verification/scheme-comparison.jl`, advisory, exits 0). The benchmark suite
  times a bare kernel and `test_convergence` measures an order on a shifted
  sine; neither says what a scheme costs *in the physics*. Ranked by error in
  the Landau damping rate at k = 0.5: `LaxWendroff` 0.64%, cubic
  `SemiLagrangian` 1.07%, `PFC` 1.45%, `Godunov`+`VanLeer` 3.99%, and upwind —
  with `Godunov(PiecewiseConstant)` and linear `SemiLagrangian`, identical to it
  as `test_amplification` requires — at **48.8%**, its own dissipation being two
  orders of magnitude larger than the physical damping it is trying to measure.

  The interesting half is the second table. A 1% perturbation cannot rank the
  schemes on positivity at all: every one returns the same `min f = 1.3e-4`,
  which is only the Maxwellian's tail at `v = ±4`. At 50% amplitude the two
  schemes that *lead* the accuracy table are exactly the two that drive `f`
  negative — `LaxWendroff` to −0.094 and cubic `SemiLagrangian` to −0.098,
  against a peak of 0.6. That is Godunov's theorem arriving in the physics, and
  it is why the harness defaults to `PFC` despite it not leading the first table.

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
  form. Six symbols are derived from the update formulas and asserted mode by
  mode: `Upwind`, `Godunov(PiecewiseConstant)` and `SemiLagrangian(LinearSpline)`
  share `1 − c(1 − e^{−iθ})`; `LaxWendroff` is `1 − ic·sinθ + c²(cosθ − 1)`;
  `Godunov(PiecewiseLinear, NoLimiter)` collapses to the centred flux
  `1 − ic·sinθ`; and `PFC`'s unlimited branch is its third-order flux. Measured
  agreement over `m ∈ {1,2,4,8,16,24,31}` and `c ∈ {0.4, 0.8}`: 8.9e-15 for the
  three-point schemes, 4.6e-14 for `PFC`.

  This is sharper than what guarded these kernels before. `test_convergence`
  fits a slope to ±0.15 and `test_golden` pins one dataset for eight steps;
  this pins every mode from the fundamental to the grid scale against an
  analytic value, and because the comparison runs elementwise it asserts at the
  same time that the output *is* a pure mode. Two consequences worth naming: the
  two claimed scheme equivalences now hold mode by mode rather than on one
  profile, and `PiecewiseLinear`'s unconditional instability is stated
  analytically — `|g| > 1` for every mode, worst 1.077 per step at `c = 0.4` —
  where it was previously demonstrated by marching `4N/c` steps and watching the
  amplitude reach 6.66e+24.

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
  match.


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
