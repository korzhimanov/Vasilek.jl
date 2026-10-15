# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project has not been released; entries below describe work on `master`.
The measurement journals the entries once carried are in the git history,
as they stood through commit `f5032f9`: `git show f5032f9:CHANGELOG.md`.

## [0.2.0] - unreleased

Each section gathers what the pull requests appended to it, newest first. An
entry says what changed and, when it needs any, where its evidence lives: a
section of `docs/src/verification.md`, `docs/src/driver-notes.md`, the migration
guide `docs/src/migration-0.2.md`, or a test.

### Breaking

- **`advect!` and `collide!` refuse a workspace narrower than the data**, as a Float32
  grid's `PFCNonUniform` or `OnGrid` and `BGK(0.5f0)` build, and `collide!` arguments that
  share memory. Migration guide, "Workspaces" and "Collisions".
- **Round-off no longer stops a run at a checked `PFC`'s bounds**: it takes a line
  within round-off of them and returns one clamped inside, so `checked` changes
  the last bits, and `Godunov` keeps a constant line (golden ≤ 4 ulps).
  `docs/src/driver-notes.md`, migration guide's "PFC".
- **`Landau1P(A)` defaults to `Tₜ = 1`, where it was 1e-3, and `Tₜ` is the
  bath's temperature**, no longer a factor on the rate; one below `Δv²` is
  refused. See the migration guide's "Collisions".
- **`Yee1D` refuses what would run wrong, and its source sees one time type**:
  bad steps, a `pml` built for other steps and a current of the wrong `axes` are
  refused; the operator copies its layer. `test/MaxwellSolver/test_fdtd_1d.jl`.
- **The FDTD solver is a type, `Yee1D`, and `make_advance_fields` is gone.**
  `advance!(mesh, op, t, j)` steps in the type of `Δx` and `Δt`, to the bit of
  the closure in `Float64`. See the migration guide's "FDTD".
- **`PoissonFourier1D.generate_solver` is gone**; `PoissonFFT1D` replaces it
  (Changed).
- **`StrangSplitting.make_time_step_2d!` is gone; `strang_step!` steps one
  matrix**, `f[x, v]`, and hands `cv` the state transposed, `f[v, x]`. See the
  migration guide's "Strang splitting".
- **`VlasovPoisson1D1V.line_advector` and `substeps` are gone**: `OnGrid` is the
  scheme on its grid, `substeps` is `Vasilek.Advection.substeps`, and a uniform
  scheme on a non-uniform grid is an `ArgumentError`.
- **Advection schemes are types, and `generate_solver` is gone.** `advect!`
  dispatches on an immutable value with an explicit `workspace`, so lines can be
  threaded. See [docs/src/migration-0.2.md](docs/src/migration-0.2.md).
- The `Symbol`-valued options are singleton types: `:Riemann_linear` becomes
  `PiecewiseLinear()`, `:VanLeer` becomes `VanLeer()`, `:Cubic` becomes
  `CubicSpline()`.
- **Collision operators follow the same shape**: `BGK(τ)` and
  `collide!(dest, src, op, v, Δt, ws)` replace `BGK.generate_solver`, which
  mutated its argument; `Landau1P` is unexported.
- `Limiters` is gone; the limiters are callable types in `Advection`. At the
  change every scheme matched 0.1 to the bit; `Godunov(PiecewiseLinear())` has
  changed since on purpose (Fixed).

### Added

- **`CLAUDE.md` and a SessionStart hook for Claude Code cloud sessions**: the
  commands and the rules for changes; `.claude/hooks/session-start.sh` installs
  the latest Julia through `juliaup` and instantiates every environment.
- **`OnGrid(scheme, Δz)`, exported: a scheme on its grid**, taking a
  displacement and splitting a step wider than a cell; `line_advector` as a
  value. `test/VlasovSolver/test_strang_splitting.jl`.
- **`strang_step!` takes a collision hook**, `Collide(op, v, Δt)`, which
  collides each v line for `Δt/2` either side of its kick.
  `test/VlasovSolver/test_strang_splitting.jl`.
- **The documentation is published, with the verification studies rendered on
  it**, at https://korzhimanov.github.io/Vasilek.jl/dev/ (`docs/make.jl`,
  `.github/workflows/Docs.yml`); every build runs the studies.
- **Collisional Landau damping: `BGK` in a run, against its own dispersion
  relation**, `PartialBGK` the control. See
  [collisional damping](docs/src/verification.md#collisional-damping).
- **The bump-on-tail instability, carried through saturation**: growth on the
  kinetic root of the grid's field, saturation by trapping, the plateau. See
  [the bump-on-tail instability](docs/src/verification.md#the-bump-on-tail-instability).
- **`Superbee`**, the most compressive limiter keeping
  `Godunov(PiecewiseLinear())` second order and TVD. See
  [the scheme comparison](docs/src/verification.md#the-scheme-comparison).
- **A nonlinear equilibrium, run to see whether it stays put**; `vlasov_poisson`
  takes the ions as `nᵢ`, and `renormalize`. See
  [the BGK equilibrium](docs/src/verification.md#the-bgk-equilibrium).
- **The plasma echo, with the field off and with it on**, against a closed form
  and against second-order kinetic theory (`test/echo.jl`). See
  [the plasma echo](docs/src/verification.md#the-plasma-echo).
- **Time reversal of the self-consistent system**, and `vlasov_poisson` returns
  `f`. See [reversibility](docs/src/verification.md#reversibility).
- **Strong Landau damping is measured against the literature**: γ₁ and γ₂ at α =
  0.5 under refinement and on the non-uniform velocity grid. See
  [strong Landau damping](docs/src/verification.md#strong-landau-damping).
- **Galilean invariance of the Landau run**: a drifting plasma damps as the one
  at rest once the Doppler phase is removed. See
  [the Doppler shift](docs/src/verification.md#the-doppler-shift).
- **The dispersion relation of an electromagnetic wave in plasma is asserted**
  (`test/test_em_plasma.jl`), with the cutoff and the current's sign; the
  current's three lines are `transverse_step!` in the harness.
- **The kinetic dispersion relation is solved rather than tabulated**:
  `Z(ζ) = i√π·erfcx(−iζ)`, `landau_root` and `two_stream_warm` in
  `test/dispersion.jl`, held to identities in `test/test_dispersion.jl`.
- **The two-stream case is measured against warm beams**, and `kv₀ = 1.0`, where
  the cold form gives zero, is its sharpest case. See
  [the two-stream instability](docs/src/verification.md#the-two-stream-instability).
- **The laser wakefield study drives a wake**, through the ponderomotive force,
  measured against `linear_wake`. See
  [the laser wakefield](docs/src/verification.md#the-laser-wakefield).
- **The two-stream instability**, the first unstable case, measured then against
  the cold closed form and since against the warm root (above), with
  `growth_rate` fitted in an amplitude band.
- **The wakefield study was asserted to run and stay bounded** (`wakefield` in
  the harness), before it had a ponderomotive coupling; the entry above replaces
  that with a measured wake.
- **The splitting and a real scheme are measured together on a rigid rotation**:
  each scheme keeps its own order. See
  [convergence](docs/src/verification.md#convergence-under-refinement).
- **The splitting is reversible, and the round trip measures dissipation**
  (`test/VlasovSolver/test_strang_splitting.jl`). See
  [reversibility](docs/src/verification.md#reversibility).
- **The Yee scheme's numerical dispersion relation**,
  `sin(ωΔt/2) = cfl·sin(kΔx/2)`, measured to 1.3e-4 from the mode's projection
  (`test/MaxwellSolver/test_fdtd_1d.jl`).
- **PML reflection against layer thickness**: four orders between 2 and 8 cells,
  the cubic ramp working (`test/MaxwellSolver/test_fdtd_1d.jl`).
- **`BGK` relaxes at the rate it is given**: `f − M` decays as `exp(−kΔt/τ)` to
  2e-10, an identity once `M` is fixed by the conserved moments.
  `test/BoltzmannSolver/test_damping_1v.jl`.
- **The schemes are compared with each other**, the L² error per scheme and
  problem class (`test/test_comparison.jl`). See
  [the scheme comparison](docs/src/verification.md#the-scheme-comparison).
- **The collision operators' complexity class is gated**: fitted exponents 1.00
  for `BGK` and 1.99 for `Landau1P`, which a uniform slowdown cannot move.
  `test/test_comparison.jl`.
- **A work–precision report** (`benchmark/workprecision.jl`, advisory): error
  against cost per scheme and problem class, with the frontier. See
  [the scheme comparison](docs/src/verification.md#the-scheme-comparison).
- **The suite runs against the Julia prerelease, weekly**, as a
  `continue-on-error` job, and every CI job runs on Mondays' cron.
  `.github/workflows/CI.yml`.
- **The verification driver takes its schemes as arguments**, `scheme_x` and
  `scheme_v`, and returns a NamedTuple. See
  [the driver notes](docs/src/driver-notes.md).
- **Mass, momentum, L² and entropy are measured and asserted**, by the
  cell-width sums the schemes conserve. See
  [invariants](docs/src/verification.md#invariants).
- **Landau damping converges under refinement**, and the L² dissipation says why
  γ sits above the root. See
  [convergence](docs/src/verification.md#convergence-under-refinement).
- **A cross-scheme study on a physical observable**
  (`verification/scheme-comparison.jl`, advisory). See
  [the scheme comparison](docs/src/verification.md#the-scheme-comparison).
- **Landau damping at three wavenumbers, and the real frequency.** See
  [Landau damping](docs/src/verification.md#landau-damping-and-dispersion).
- **The frequency estimator is held to its window**, as `γ` is. See
  [Landau damping](docs/src/verification.md#landau-damping-and-dispersion).
- **The plasma oscillation frequency is asserted against Bohm–Gross**, and the
  cold `ωₚ` is excluded. See
  [plasma oscillations](docs/src/verification.md#plasma-oscillation-frequency-and-energy).
- **The extended verification suite runs in CI**, on every pull request, in the
  `Extended verification` job of `.github/workflows/CI.yml`.
- **Von Neumann amplification factors, against closed forms**, mode by mode for
  six schemes, to 1.8e-14 (`test/test_amplification.jl`).
- **Dissipation and dispersion, per scheme per wavenumber**, from the same
  symbols (`test/test_amplification.jl`).
- **Free streaming: phase mixing and recurrence**
  (`test/VlasovSolver/test_free_streaming.jl`). See
  [recurrence](docs/src/verification.md#recurrence).
- **The Yee leapfrog's conserved energy is asserted**, the staggered
  `‖E^{n+1}‖² + ⟨H^{n+1/2}, H^{n+3/2}⟩`, to 3e-15 over 2000 steps; the naive
  `‖E‖² + ‖H‖²` is not conserved. `test/MaxwellSolver/test_fdtd_1d.jl`.
- **PML absorption is measured**: reflection 6.25e-9 either way, against a
  control without the layer. `test/MaxwellSolver/test_fdtd_1d.jl`.
- **`BGK` satisfies the H-theorem**, to machine precision once the velocity
  window resolves the relaxed state. `test/BoltzmannSolver/test_damping_1v.jl`.
- **`∂f∂v` is tested directly**: second order inside, first at the ends, exact
  on a linear profile, on a non-uniform grid too.
  `test/BoltzmannSolver/test_damping_1v.jl`.
- `YeeMesh1D` shape and element type, and a `Bounds checking` CI job, since
  replaced by `@inbounds in effect` (Changed).
- **The README's usage example is executed by the suite**
  (`test/test_readme.jl`), and so cannot drift.
- **Reentrancy is asserted**: lines swept by tasks sharing a scheme value
  reproduce the serial result bit for bit, run by a `Multithreaded` CI job.
  `test/test_threading.jl`.
- **`BGK` moment conservation**, on skewed and drifting data, an arbitrary
  Maxwellian as a fixed point, and both limits of the update.
  `test/BoltzmannSolver/test_damping_1v.jl`.
- **Contract tests** for `advect!`: argument checks, stateless workspaces,
  `@inferred` on every kernel, element types. `test/test_contracts.jl`.
- **Direction-symmetry suite**: `R∘A(+c) == A(−c)∘R` for every scheme, to
  4.5e-14 at worst, and `SemiLagrangian`'s periodicity in `c`.
  `test/test_symmetry.jl`.
- **`StrangSplitting` has tests**: second order in Δt on a rigid rotation, and
  the transpose bookkeeping. `test/VlasovSolver/test_strang_splitting.jl`.
- **`PFC`'s bounds check is tested**, at both ends and both branches, and
  `checked = false` is bit for bit the same on valid data.
  `test/test_contracts.jl`.
- **The Poisson solver's exact discrete identity is asserted**,
  `E_num == sinc(kΔx)·E_exact` to 1.2e-13.
  `test/MaxwellSolver/test_poisson_fourier_1d.jl`.
- **FDTD at the magic time step** translates a pulse exactly, to 1.1e-19, and
  the Courant limit is asserted. `test/MaxwellSolver/test_fdtd_1d.jl`.
- **Extended verification** (`test/test_verification.jl`), behind
  `VASILEK_EXTENDED=1`: the studies' claims as assertions. See
  [what the tests assert](docs/src/verification.md#what-the-tests-assert).
- [docs/src/normalization.md](docs/src/normalization.md): the units, the
  wavenumber and FDTD current conventions, and why PFC's bounds have no default.
- `PFC(; fmin, fmax, checked = false)` compiles the bounds check away; the check
  is a type parameter, on by default.
- **Allocation gate**, blocking and deterministic, at zero for every kernel that
  achieves it. `test/test_allocations.jl`.
- `--strict` and `--rebaseline` flags for the benchmark driver, and
  `PFCNonUniform` and `FDTD1D` in the benchmark suite.
- **Golden-value regression tests**, bit for bit, for every advection scheme;
  regenerate with `test/generate_golden.jl` only for an intended change.
  `test/test_golden.jl`.
- **Order-of-accuracy suite**, each scheme's global order measured before being
  asserted. `test/test_convergence.jl`.
- **Structural invariants**: mass, total variation, PFC positivity and the
  maximum principle, and the Courant limit. `test/test_invariants.jl`.
- Continuous integration on GitHub Actions: Julia `lts` and `1` on Linux,
  Windows and macOS, and a coverage job reporting to Codecov.
- Dependabot updates for GitHub Actions.
- `Aqua.test_all` in the test suite.
- Separate `benchmark/Project.toml` and `verification/Project.toml`.
- Advection tests asserted that both `generate_solver` overloads gave the same
  bits, the invariant `PFC` broke; the overloads have since gone (Breaking).

### Changed

- **The API reference is a page per module, and a docs page over 100 KiB fails
  the build**: the one-page API had reached Documenter's warning limit. Old
  `api/#Name` anchors now land on the overview. `docs/make.jl`.
- **The normalization notes link three measurements to the verification page**,
  and the README links the built notes; the harness's `em_omega` and `wakefield`
  docstrings keep the carrier's 4.5% apart from the pulse's 5.2% and 6.6%.
- **The measurements have one home**: `docs/src/verification.md` gathers them by
  study; the README keeps ten rows of its table, held to it by
  `test/test_readme.jl`, and this file one to three lines an entry.
- **`Landau1P` is a consistent model**: a regularised kernel, the transverse
  bath's drag, and a flux-form update that keeps the density to 2.4e-16; its
  Known issues are closed. `test/BoltzmannSolver/test_damping_1v.jl`.
- **The quadratic and cubic `SemiLagrangian` prefilter in the package**, ten
  times faster, allocating nothing, and moving the splines by 4.4e-15 at most.
  `test/VlasovSolver/test_spline.jl`.
- **The two polarisations of the FDTD step share one update**, identical in
  every byte and two to five times faster. `test/MaxwellSolver/test_fdtd_1d.jl`.
- **A step of `vlasov_poisson` allocates nothing**, every Float64 result the
  same to the bit, and a run continues by calling the driver on the `f` it
  returned. `test/test_allocations.jl`, `test/test_driver.jl`.
- **`vlasov_poisson` steps through `strang_step!`, to the bit**, every line
  contiguous; under Float32 a line stays in Float64 through its sub-steps and
  its kick. `test/test_driver.jl`.
- **Uniform-grid schemes are split into sub-steps too**, by `OnGrid`: a step
  within a cell is unchanged to the bit, and one past `|c| = 1` now runs.
- **`test/VlasovSolver/test_strang_splitting.jl` steps through `strang_step!`**,
  holding the step, the hook and the sub-steps to their compositions, bit for
  bit.
- **The docs build with Documenter**, from `docs/src/`, with an index and a page
  of every docstring (`docs/make.jl`).
- **`Pkg.test(test_args = ["golden", "contracts"])` runs those files only**;
  `test/runtests.jl` lists every file and refuses an unknown name.
- **`vlasov_poisson` is held to its step bit for bit**, the step written out
  from `advect!`, `collide!`, `PoissonFFT1D` and the cell widths.
  `test/test_driver.jl`.
- The centred Poisson difference and `slope_limit` lose their Float64 literals,
  the same bits in Float64.
- `Scripts.yml` runs on a push to master only for the paths its pull requests
  watch.
- `runbenchmarks.jl` prints the baseline's Julia and BenchmarkTools versions
  beside the run's.
- The complexity-class test samples for longer, and a test holds the manual's
  front page to the README's introduction (`test/test_docs.jl`).
- **Documentation says what the code does**: docstrings state contracts, the
  driver's reasons are in `docs/src/driver-notes.md`, and `test/test_docs.jl`
  runs the code in `docs/src/*.md`. The version is `0.2.0-DEV` until a release.
- **CI runs what it says it runs**: an `@inbounds in effect` job, the extended
  verification on the LTS too and failing without its marker file, coverage in
  one job, declared permissions. `.github/workflows/CI.yml`.
- **The scripts are run**: `.github/workflows/Scripts.yml` runs every study to
  completion and loads the benchmark suites.
- **CompatHelper watches the `[compat]` bounds** of the package and of the
  script environments.
- **The script environments install from a clean checkout**, through
  `[sources]`, which puts their `julia` floor at 1.11.
- **Benchmarks** use `N` points rather than `N + 1`, which timed FFTW's slow
  path at prime lengths, and report entries untuned or without a baseline.
- **`BGK` relaxes to the discrete Maxwellian and conserves to round-off** on any
  window; the sampled one is `BGK(τ; conservative = false)`. See
  [collisional damping](docs/src/verification.md#collisional-damping).
- **`NumericalIntegration` is no longer a dependency**, only a test one.
- **The 1D1V driver ships with the package: `vlasov_poisson`, exported**, in
  `Vasilek.VlasovPoisson1D1V`; it refuses a non-uniform x grid and fewer than
  two times. See [the driver notes](docs/src/driver-notes.md).
- **`verification/landau-damping-1d1v.jl` and `plasma-oscillations-1d1v.jl` run
  the driver**, and take their spectra from its mode amplitude, where they kept
  five copies of the loop.
- **The README says what exists**, 1D1V on static grids, and shows a whole
  driver run.
- **The Poisson solver is a value, `PoissonFFT1D(n, Δx)`, with a `workspace`**;
  `:centered` is 0.1's bit for bit, `:spectral` exact for every resolved mode.
  See the migration guide's "Poisson".
- **`strang_step!` took schemes and a workspace**, the closure form's step to
  the bit; it has since become the split step itself (Breaking).
- **`make_time_step_2d!` left both arrays holding the step**; it has since gone
  (Breaking).
- **`workspace` is defined once, in `Vasilek`**, and every module adds methods
  to it.
- **`PFC`'s and `PFCNonUniform`'s bounds check throws `DomainError`**, not
  `AssertionError`, which may be compiled out.
- **`advect!` refuses more wrong calls**: a `dest` that may alias `src`, offset
  arrays, a non-finite Courant number. `test/test_contracts.jl`.
- **Constructors refuse what the schemes cannot use**: `fmin > fmax`, bad cell
  widths for `PFCNonUniform`, a limiter with `PiecewiseConstant()`.
  `test/test_contracts.jl`.
- **Element types are carried through**: `PFCNonUniform` on Float32 data,
  `workspace(scheme, n, T)`, no Float64 literals in the kernels; Float64 results
  are bit-identical.
- **`Landau1P` is no longer exported from `Collisions`**; reach it as
  `Vasilek.Collisions.Landau1P`.
- **`Godunov(PiecewiseLinear(), VanLeer())` steps 7 to 14 times faster and
  `Superbee` 6.5 to 11.5 times, to the same bits**: the flux inlined, `iszero`
  for `≈ 0.0`, each face evaluated once. `test/test_golden.jl`, `benchmark/`.
- **`FDTD1D` drives only the interior nodes, and the PEC boundaries are written
  down**: a pulse reaching either wall returns inverted, at −0.9998.
  `test/MaxwellSolver/test_fdtd_1d.jl`, `docs/src/normalization.md`.
- **`advect!` validates its arguments**: an aliased `dest`, unequal lengths, a
  workspace of the wrong size, fewer than three cells. `test/test_contracts.jl`.
- **The verification notebooks are runnable scripts**, Literate `.jl` in place
  of `.jmd` and Weave.
- **`wakefield.jl` produces output and is stable**, at a peak field of 0.38
  where the current's bugs gave 1.0e22. See
  [the laser wakefield](docs/src/verification.md#the-laser-wakefield).
- `Upwind` and `LaxWendroff` reduced to a single numerical kernel each, the same
  bits, confirmed by the golden tests.
- Tests load the package (`using Vasilek`) instead of including source files by
  relative path.
- Runtime dependencies cut from nine to four, and `[compat]`, `[extras]` and
  `[targets]` added.
- `@fastmath` removed from the advection kernels: it bought nothing, and it
  licensed the reassociation that bit-for-bit tests rule out.
- `Limiters` hoisted to the top level, so the test process holds one copy.
- `StrangSplitting` added to the package include list.
- `Manifest.toml` is gitignored rather than tracked, as for a library.

### Removed

- **Interpolations as a dependency of the package**; it is a test dependency,
  for `test/VlasovSolver/test_spline.jl`. The package depends on `FFTW` and
  `LinearAlgebra`.
- `StrangSplitting.make_time_step_2d!`, `VlasovPoisson1D1V.line_advector` and
  `VlasovPoisson1D1V.substeps` (now `Advection.substeps`); see Breaking.
- Five files holding a no-op function each and the three untyped `src/Mesh`
  types, referenced by nothing.
- The `plot_needed` / `plot!` machinery in tests, the only reason the test
  environment needed Plots.

### Fixed

- **`OnGrid` takes a grid uniform to rounding, as `make_poisson` does**, with `scale`
  for widths taken from coordinates; positional `PFC` constructors round bounds outward.
  `test/test_driver.jl`; migration guide, "Advection" and "Note on the fourth argument".
- **`Godunov(PiecewiseLinear(), VanLeer())` stays finite where the ratio of two
  finite differences overflows**: `VanLeer` is 2 at `r = Inf`, as `AbstractLimiter`
  now asks of every limiter, and else the same to the bit. `test/test_contracts.jl`.
- **`vlasov_poisson` runs in the data's type, `T = float(eltype(f₀))`**: field,
  density, histories and `E_modes` follow it, and the lines are worked in
  Float64. `test/test_driver.jl`.
- **`PFC` keeps a line that sits at its upper bound**: a cell takes
  `f + (Φin − Φout)`, zero on a constant line, so the next checked call no
  longer refuses it. `test/test_contracts.jl`.
- **The default schemes keep a line that sits at their upper bound**:
  `PFCNonUniform` takes `f + (Φin − Φout)/Δx`, and a waterbag stays an
  equilibrium to the bit. `test/VlasovSolver/test_nonuniform_advection.jl`.
- **`cell_widths` keeps a Float32 grid Float32.**
- **`vlasov_poisson` runs grids and data of other types**: the widths follow the
  data, the default schemes work in Float64 at least, and Float64 runs are as
  they were to the bit. `test/test_driver.jl`.
- **`PFCNonUniform` rounds its bounds outward** to its element type, and
  `keeps_bounds` compares in it.
- **`PFC{T, Checked}` and `PFCNonUniform{T, Checked}` refuse a `Checked` that is
  not a Bool**, and check what the keyword constructors check.
- **`collide!` refuses a `src` sharing its workspace**, which returned the
  Maxwellian alone.
- **A scheme that takes no workspace refuses one.**
- **`advect!` refuses a workspace whose buffer shares memory with `dest` or
  `src`.** `test/test_contracts.jl`.
- **Tests that could not fail now test the code**: the `Landau1P` antisymmetry,
  scheme values against their `deepcopy`, a root against itself, the dissipation
  symbols against their own closed forms.
- **The extended suite's skip is a `@test_skip`**, so a run that lost
  `VASILEK_EXTENDED` shows Broken.
- **Fragile checks loosened to what they mean**: the allocation gate allows one
  boxed value between N and 4N. The spline golden rows, loosened then, are held
  to `==` since (Changed).
- **The harness**: `zero_crossings` counted an exact zero twice, `wakefield`
  drove the laser one Δt ahead of its clock, and the reversibility test
  normalised by its own result.
- **The driver's docstring says when the invariants are sampled**: at the end of
  each step, half a step after the field.
- The test files calling `norm` import it, so each runs on its own.
- **The z-polarised FDTD source launched its pulse into the absorbing layer**;
  both polarisations launch rightwards. `test/MaxwellSolver/test_fdtd_1d.jl`.
- **`workspace(BGK(τ), n)` was a `MethodError`**; the collision operators extend
  `Vasilek.workspace`.
- **`BGK` filled an empty or cold velocity line with NaN**; a line colder than
  the grid resolves, `T < Δv²`, passes through unchanged.
- **`PML` coefficients**: `r₂` is taken through `expm1`, and a negative `N` or
  `σ_max` is refused.
- **`make_advance_fields` refused inconsistent arguments**; `Yee1D` has since
  replaced it (Breaking).
- **The harness's Poisson solve threw on about a quarter of grids**; it takes
  `rfftfreq`, tested over every `Nx` in 8:512.
- **The verification environment can load the harness again**, with
  `NumericalIntegration` back in `verification/Project.toml`.
- **The charge rescaling is the sums the solver conserves**, `Σ f ΔvΔx` and
  `Σ nᵢ Δx`: 1 to round-off, where the trapezoid scaled `ωₚ²` by up to 7.9e-3.
  See [the driver notes](docs/src/driver-notes.md).
- **What the charge fix moved, and quotes stale before it**, corrected
  throughout, and a test pins the strong-damping seed's mass to 4e-15; the
  figures before and after are in the git history.
- **The energy histories are the sums the schemes conserve, at one instant**:
  the kinetic energy is centred on the kick, so `ε` is sampled at `t[k] + Δt/2`
  as `ε_e` is. See [invariants](docs/src/verification.md#invariants).
- **`PFCNonUniform` checks its data against its bounds, as `PFC` does**, for 5%
  of a step on Julia 1.13; data outside them was corrupted without an error.
  `test/test_contracts.jl`.
- **A non-positive scheme given for one direction, the other the default, is
  refused before the first step**: `keeps_bounds` decides. See
  [the driver notes](docs/src/driver-notes.md).
- **With `collisions`, the driver's defaults lift their upper bound**, since BGK
  raises a line past its maximum. See
  [collisional damping](docs/src/verification.md#collisional-damping).
- **`two_stream(0.4)` ran on 95 cells rather than 96**, a float range coming up
  a point short; periodic grids are built by length. See
  [the two-stream instability](docs/src/verification.md#the-two-stream-instability).
- **`Godunov(PiecewiseLinear())` is second order, and TVD to `|c| = 1`**: its
  flux gained the strip average's `(1 − |c|)`, which also makes the unlimited
  scheme `LaxWendroff`. `test/test_convergence.jl`, `test/test_invariants.jl`.
- **`advect!` refuses a step past the scheme's Courant limit**, a `DomainError`
  before anything is written; `SemiLagrangian` opts out, and `PFCNonUniform` is
  held to its narrowest cell. `test/test_invariants.jl`.
- **A full test run evaluates each shared test file once**: each include is
  guarded on a name the file defines, and a full run prints none of the 150
  method-overwrite warnings.
- **`test_free_streaming.jl` replaced the harness's `mode_amplitude`** for every
  file run after it; its own helper is `cos_amplitude`.
- **`PFC` in the remaining verification runs is bounded by the distribution it
  carries**: the driver, `ballistic_echo` and `free_stream` take a scheme as a
  function of the initial `f`.
- **`docs/src/normalization.md` gave the Poisson equation with the wrong sign**;
  it is `∂E/∂x = nₑ − nᵢ`, as every solve takes it. See
  [the BGK equilibrium](docs/src/verification.md#the-bgk-equilibrium).
- **The verification harness gave PFC an upper bound of 1 whatever `f` was**;
  the bounds are those of `f` on entry, and every bounded run stays inside them.
  See [the BGK equilibrium](docs/src/verification.md#the-bgk-equilibrium).
- **The "linear" Landau case at k = 0.3 was not linear**, and its agreement was
  two errors cancelling; the cases run at α = 1e-3, guarded by `trapping_phase`.
  See [Landau damping](docs/src/verification.md#landau-damping-and-dispersion).
- **The recurrence in the Landau notebook is the second harmonic's**, and
  `vlasov_poisson` records modes through `modes`. See
  [recurrence](docs/src/verification.md#recurrence).
- **The laser in the wakefield study was travelling six percent slow**, the Yee
  grid's dispersion, which `vg_pulse` predicts. See
  [the laser wakefield](docs/src/verification.md#the-laser-wakefield).
- **"Finite temperature cannot make a beam grow faster than a cold one" is
  false**: `a = 0.8` overshoots by physics, the `vt` sweep only under a fixed
  window. See [two-stream](docs/src/verification.md#the-two-stream-instability).
- **The Landau damping rate was fitted with an estimator its own window chose
  the answer for**; it is fitted through the maxima. See
  [Landau damping](docs/src/verification.md#landau-damping-and-dispersion).
- **`local_extrema` skipped the last interior point it was asked about.**
  `test/test_harness.jl`.
- **The `(h, h₀)` methods of `LaxWendroff` sized their loops from the captured
  array**, as `PFC`'s did; the closures have gone since (Breaking).
- `Godunov.generate_Φ` and `Limiters.generate_limiter` returned `nothing` for an
  unrecognised option; both have gone since (Breaking).
- **`SemiLagrangian` rebuilt its interpolation on every step**; it prefiltered a
  buffer in place, and the package's own prefilter has since taken the splines
  to zero allocations (Changed).
- **`BGK` allocated six temporaries per step**; it allocates nothing, to the
  same bits.
- Benchmark measurements were invalid: `@benchmarkable` timed a global `Dict`
  lookup. The closures are interpolated.
- The benchmark driver ran only the Vlasov group; `results.json` stores the
  estimate, 13 kB, where it stored 4.4 MB of samples.
- The Cyrillic `c` (U+0441) in the benchmark keys is gone.
- **`FDTD1D` read `f.hz[i-1]` in an `hy`-only stencil for `ez`**, so the y
  polarization drove `ez`: numerically breaking for any z-polarized run
  (`test/MaxwellSolver/test_fdtd_1d.jl`).
- `FDTD1D`'s `PML` default argument passed `Δt` and `Δx` in the wrong order; a
  keyword constructor makes that unrepresentable.
- **`Landau1P` differentiated at the wrong index**, and integrated an unassigned
  `J[i]`; `L` and `Tₜ` are keywords (`test/BoltzmannSolver/test_damping_1v.jl`,
  "Landau1P invariants").
- **`PFCNonUniform` computed its limiter from the global spacing ratio**; per
  cell, the non-uniform plasma oscillation drifts 4.75%, not 12%. See
  [plasma oscillations](docs/src/verification.md#plasma-oscillation-frequency-and-energy).
- `PFCNonUniform` hardcoded `fₘᵢₙ = 0`, `fₘₐₓ = 1`; they are required keywords,
  as in `PFC`. A breaking change.
- **`PoissonFourier1D` returned a field too large by `1/Δx²`**: numerically
  breaking for anything computed with it before
  (`test/MaxwellSolver/test_poisson_fourier_1d.jl`).
- **`PFC`'s two `generate_solver` overloads computed different things**, a
  factor of 740 on data above 1; `fₘᵢₙ` and `fₘₐₓ` are required. A breaking
  change.
- `PFC`'s three-argument `solve!(h, h₀, c)` sized its loops from the captured
  array.
- The Poisson test sampled 1001 points of a period-10 sine, which leaked
  spectrally.
- `test_damping_1v.jl` rebound a local in its stepping loop, so its "100 time
  steps" repeated step one; the `Landau1P` bug the running loop exposed is fixed
  (Changed).
- `test_fdtd_1d_current` had its assertion commented out; it is restored and
  passes.
