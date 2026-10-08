#!/bin/bash
# Installs the latest stable Julia through juliaup and instantiates the
# package's environments, so that tests and verification studies run at once
# in a Claude Code cloud session. Idempotent; does nothing outside the cloud.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

export PATH="$HOME/.juliaup/bin:$PATH"

# `release` is the current stable; `update` moves it to a newer one if there is.
# The steps are chained with `&&` because `set -e` does not apply inside a
# function called as an `if` condition. Their output goes to stderr: a
# SessionStart hook's stdout lands in the session's context.
install_julia() {
  { command -v juliaup >/dev/null 2>&1 ||
    curl -fsSL https://install.julialang.org | sh -s -- --yes --add-to-path=no; } &&
  juliaup add release &&
  juliaup update release &&
  juliaup default release
}
# A failed update leaves an installed Julia usable, so only its absence stops
# the hook; either way the session starts.
if ! install_julia >&2; then
  echo "session-start: Julia could not be installed or updated (are" \
       "install.julialang.org and julialang-s3.julialang.org allowed by the" \
       "network policy?)" >&2
  command -v julia >/dev/null 2>&1 || exit 0
fi

echo "export PATH=\"$HOME/.juliaup/bin:\$PATH\"" >> "${CLAUDE_ENV_FILE:-/dev/null}"

# The repository root, also when the hook is run by hand to debug it.
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}"

# The package and its test-only dependencies ([extras] are not in the
# environment itself, so they are added to the depot by name), then the
# verification studies (Plots included) and the documentation. Packages come
# from pkg.julialang.org, which redirects to its regional mirrors and to
# storage.julialang.net; if the network policy blocks them, say so and let the
# session start with Julia alone rather than fail it.
instantiate() {
  julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' &&
  julia -e 'using Pkg; Pkg.activate(; temp = true);
            Pkg.add(["Aqua", "NumericalIntegration", "SpecialFunctions"])' &&
  julia --project=verification -e 'using Pkg; Pkg.instantiate()' &&
  julia --project=docs -e 'using Pkg; Pkg.instantiate()'
}
if ! instantiate >/dev/null 2>&1; then
  echo "session-start: Julia is installed, but the packages could not be" \
       "instantiated (are pkg.julialang.org, *.pkg.julialang.org and storage.julialang.net allowed by the network policy?)" >&2
fi
