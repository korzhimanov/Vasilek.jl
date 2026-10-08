#!/bin/bash
# Installs the latest stable Julia through juliaup and instantiates the
# package's environments, so that tests and verification studies run at once
# in a Claude Code cloud session. Idempotent; does nothing outside the cloud.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

export PATH="$HOME/.juliaup/bin:$PATH"

# Later Bash commands of the session find juliaup's Julia; the line is written
# once, however many times the hook runs.
path_line="export PATH=\"$HOME/.juliaup/bin:\$PATH\""
if [ -n "${CLAUDE_ENV_FILE:-}" ] &&
   ! grep -qxF "$path_line" "$CLAUDE_ENV_FILE" 2>/dev/null; then
  echo "$path_line" >> "$CLAUDE_ENV_FILE"
fi

# The hook's input names what started the session. After /clear or a
# compaction the container is the one set up at startup or resume, so nothing
# is installed again. Run by hand from a terminal, the hook has no input.
source=""
if [ ! -t 0 ]; then
  source=$(sed -n 's/.*"source" *: *"\([a-z]*\)".*/\1/p' || true)
fi
if [ "$source" = "clear" ] || [ "$source" = "compact" ]; then
  exit 0
fi

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
# A failed update leaves an installed Julia usable, so only its absence (or a
# Julia below the package's 1.10 floor) stops the hook; either way the session
# starts.
if ! install_julia >&2; then
  echo "session-start: Julia could not be installed or updated (are" \
       "install.julialang.org and julialang-s3.julialang.org allowed by the" \
       "network policy?)" >&2
  if ! julia -e 'exit(VERSION >= v"1.10" ? 0 : 1)' >/dev/null 2>&1; then
    echo "session-start: no Julia 1.10 or newer to instantiate the" \
         "package with" >&2
    exit 0
  fi
fi

# The repository root, also when the hook is run by hand to debug it.
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}"

# The package, then its test environment, then the verification studies
# (Plots included) and the documentation. The test environment is the
# package's [deps] and its [targets].test, under the package's [compat], so
# that `Pkg.test()` resolves the versions put in the depot here. Packages come
# from pkg.julialang.org, which redirects to its regional mirrors and to
# storage.julialang.net; if the network policy blocks them, say so and let the
# session start with Julia alone rather than fail it.
instantiate() {
  julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' &&
  julia - <<'EOF' &&
using Pkg, TOML
project = TOML.parsefile("Project.toml")
uuids = merge(get(project, "deps", Dict()), get(project, "extras", Dict()))
names = union(keys(get(project, "deps", Dict())), project["targets"]["test"])
compat = filter(kv -> kv[1] in names || kv[1] == "julia", get(project, "compat", Dict()))
env = mktempdir()
open(joinpath(env, "Project.toml"), "w") do io
    TOML.print(io, Dict("deps" => Dict(n => uuids[n] for n in names), "compat" => compat))
end
Pkg.activate(env)
Pkg.instantiate()
EOF
  julia --project=verification -e 'using Pkg; Pkg.instantiate()' &&
  julia --project=docs -e 'using Pkg; Pkg.instantiate()'
}
if ! instantiate >&2; then
  echo "session-start: Julia is installed, but the packages could not be" \
       "instantiated; the error is above (if it is a download, are" \
       "pkg.julialang.org, *.pkg.julialang.org and storage.julialang.net" \
       "allowed by the network policy?)" >&2
fi
