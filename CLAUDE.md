# CLAUDE.md

## Checking changes

- Check by actually running the code wherever possible (tests, examples, scripts
  from `verification/`, benchmarks), not only by reading the sources.
- Run with the latest stable Julia. If Julia is missing or outdated, install the
  current one (for instance through `juliaup`:
  `curl -fsSL https://install.julialang.org | sh -s -- --yes`, then
  `juliaup add release && juliaup default release`). In cloud sessions
  `.claude/hooks/session-start.sh` does this.

## Commands

- Quick tests: `julia --project=. -e 'using Pkg; Pkg.test()'`.
  Single files: `Pkg.test(test_args = ["golden", "contracts"])` — names without
  `test_` and `.jl`, listed in `test/runtests.jl`.
- Extended verification (as CI runs on a PR):
  `VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'`.
- Studies: once `julia --project=verification -e 'using Pkg; Pkg.instantiate()'`,
  then `julia --project=verification verification/<study>.jl` (needs Julia ≥ 1.11).
- Documentation: `julia --project=docs docs/make.jl`. The examples in the README
  and docs are checked by tests (`test_readme.jl`, `test_docs.jl`).

## Compatibility

- The floor is Julia 1.10 (LTS); CI runs `lts` and `1` on Linux, Windows and macOS.
  Work on the latest version, but use no language or Pkg features missing from
  1.10 (except in the `verification/` and `docs/` environments, whose floor is 1.11).
- The Manifest is not committed.

## Code requirements

- A solver step allocates nothing (`test/test_allocations.jl`) and works in the
  data's type, not only in `Float64`.
- Every accuracy claim from `verification/` is pinned by a test; the tolerance and
  the measured value go into the table in `docs/src/verification.md`.
- Golden data are never edited by hand — they are regenerated with
  `test/generate_golden.jl`.
- Units and normalization — `docs/src/normalization.md`.

## Conventions

- Every PR adds an entry to `CHANGELOG.md` (Keep a Changelog, the current
  version's section) with the measurements that justify the change. Breaking
  changes are also described in `docs/src/migration-0.2.md`.
- Commits: `fix:` / `refactor:` and so on, worded as behavior
  ("a step of vlasov_poisson allocates nothing").
- Bring a branch up to date with master by merge, not by rebase.
- Code, comments, documentation and commits are in English; talk to the user
  in Russian.
