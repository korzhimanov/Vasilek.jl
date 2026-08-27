# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project has not been released; entries below describe work on `master`.

## [Unreleased]

### Added

- Continuous integration on GitHub Actions: tests on Julia `lts` and `1` across
  Linux, Windows and macOS, plus a coverage job reporting to Codecov.
  Verified locally on both ends of the matrix (1.10.12 and 1.12.7).
- Dependabot updates for GitHub Actions.
- `Aqua.test_all` in the test suite.
- Separate `benchmark/Project.toml` and `verification/Project.toml`.
- Advection tests now assert that both `generate_solver` overloads produce
  bit-identical output — the invariant `PFC` had been violating.

### Fixed

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

- Tests load the package (`using Vasilek`) instead of including source files by
  relative path and dispatching through `eval`. `src/Vasilek.jl` was previously
  never exercised by its own test suite.
- Runtime dependencies cut from nine to four: `FFTW`, `Interpolations`,
  `LinearAlgebra`, `NumericalIntegration`. `DSP` and `Weave` were imported
  nowhere; `Plots` only inside dead branches; `BenchmarkTools` belongs to
  benchmarks. Added the missing `[compat]` section and `[extras]`/`[targets]`.
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

- `Landau1P` computes `Δfⱼ` but branches on `i` and indexes `f₀[i±1]`, so
  `Δfⱼ ≡ Δfᵢ`; and `J[j]` is never assigned when `i == j`. Marked `@test_broken`.
- `FDTD1D` uses `f.hz[i-1]` in an `hy`-only stencil for `ez`, injecting spurious
  `ez` on every run. No test covers the z polarization.
- `FDTD1D`'s `PML` default argument passes `Δt` and `Δx` in the wrong order.
  Dead code today: every call site passes a PML explicitly.
- Benchmark timings are invalid: `@benchmarkable` interpolates only `N`, so a
  non-`const` global `Dict` is looked up inside the timed region.
