#!/bin/bash
# Installs the latest stable Julia through juliaup and instantiates the
# package's environments, so that tests and verification studies run at once
# in a Claude Code cloud session. Idempotent; does nothing outside the cloud.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

export PATH="$HOME/.juliaup/bin:$PATH"

if ! command -v juliaup >/dev/null 2>&1; then
  curl -fsSL https://install.julialang.org | sh -s -- --yes --add-to-path=no
fi

# `release` is the current stable; `update` moves it to a newer one if there is.
juliaup add release
juliaup update release
juliaup default release

echo "export PATH=\"$HOME/.juliaup/bin:\$PATH\"" >> "${CLAUDE_ENV_FILE:-/dev/null}"

cd "$CLAUDE_PROJECT_DIR"

# The package and its test-only dependencies ([extras] are not in the
# environment itself, so they are added to the depot by name), then the
# verification studies (Plots included) and the documentation. Packages come
# from pkg.julialang.org; if the network policy blocks it, say so and let the
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
       "instantiated (is pkg.julialang.org allowed by the network policy?)" >&2
fi
