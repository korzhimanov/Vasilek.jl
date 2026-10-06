# Builds the documentation into docs/build: the manual in docs/src, the API, and
# the verification studies, rendered by Literate.jl with their output and
# figures. From the repository root:
#
#     julia --project=docs -e 'using Pkg; Pkg.instantiate()'           # once, Julia ≥ 1.11
#     julia --project=verification -e 'using Pkg; Pkg.instantiate()'   # once
#     julia --project=docs docs/make.jl                                # studies not run
#     VASILEK_DOCS_EXECUTE=1 julia --project=docs docs/make.jl         # studies run
#
# Running the studies is what puts their numbers and figures on the pages; it
# takes a few minutes, so a local build skips it unless asked, as the extended
# test suite does without VASILEK_EXTENDED. CI always runs them (Docs.yml), on
# pull requests too: a study can break under Literate while running fine as a
# script. `deploydocs` publishes only a build that ran them.

ENV["GKSwstype"] = "100"     # GR without a display, set before Plots loads it

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC = joinpath(@__DIR__, "src")
const EXECUTE = get(ENV, "VASILEK_DOCS_EXECUTE", "0") == "1"

# The studies load Plots and the rest from the verification environment, stacked
# under this one, so that a package a study gains there is the one it gets here.
insert!(LOAD_PATH, 2, joinpath(ROOT, "verification"))

using Documenter, Literate
using Vasilek

# ------------------------------------------------------------------- studies
# The order of the navigation: the linear electrostatic cases, the
# instabilities, the kinetic and nonlinear ones, collisions, the electromagnetic
# one, and the advisory comparison. Each page takes its title from the script's
# `# # Title` line.
const STUDIES = ["landau-damping-1d1v", "plasma-oscillations-1d1v", "two-stream",
                 "bump-on-tail", "plasma-echo", "bgk-equilibrium",
                 "collisional-damping", "wakefield", "scheme-comparison"]

# A study missing from the list would never be rendered, and nothing else would
# notice.
let scripts = sort([first(splitext(f)) for f in readdir(joinpath(ROOT, "verification"))
                    if endswith(f, ".jl")])
    sort(STUDIES) == scripts || error("docs/make.jl lists the studies $(sort(STUDIES)), " *
                                      "but verification/ holds $scripts")
end

const STUDIES_DIR = joinpath(SRC, "studies")
rm(STUDIES_DIR; recursive = true, force = true)   # nothing left over from an earlier build

# A study built without running it has no figures, and Documenter refuses an
# image whose file does not exist: the `#md # ![...](...)` lines go, and a note
# under the title says why the page is bare.
const UNEXECUTED_NOTE = """
#
# !!! note "Built without running the study"
#     This page shows the study's code but not its output or figures. The
#     published site runs every study; set `VASILEK_DOCS_EXECUTE=1` to build
#     it that way locally.
"""
function unexecuted(text)
    text = replace(text, r"^#md # !\[.*\n"m => "")
    m = match(r"^# # .*\n"m, text)
    m === nothing && error("a study has no `# # Title` line to put the note under")
    i = m.offset + ncodeunits(m.match)
    return text[1:i-1] * UNEXECUTED_NOTE * text[i:end]
end

# What Literate writes that says nothing on a page. It shows a chunk's value,
# and a chunk that ends by defining a function has one:
# `f (generic function with 1 method)`. And unexecuted, it appends
# `nothing #hide` to a chunk ending in `;`, for Documenter to run; nothing does.
tidy(text) = replace(text,
    r"\n````\n\S+ \(generic function with \d+ methods?\)\n````\n" => "\n",
    r"^nothing #hide\n"m => "")

for stem in STUDIES
    Literate.markdown(joinpath(ROOT, "verification", "$stem.jl"), STUDIES_DIR;
        flavor = Literate.DocumenterFlavor(),   # which also writes the EditURL
        execute = EXECUTE,
        # A plain fence either way: unexecuted, Literate's default for
        # Documenter is ```@example, which Documenter would then run.
        codefence = "````julia" => "````",
        # Unset, Literate asks the remote with `git remote show origin`, for
        # every page, which branch to name in links this site does not use.
        edit_commit = "master",
        preprocess = EXECUTE ? identity : unexecuted,
        postprocess = tidy)
end

# ------------------------------------------------------------------- the site
makedocs(;
    sitename = "Vasilek.jl",
    authors = "Artem Korzhimanov",
    modules = [Vasilek, Vasilek.Advection, Vasilek.StrangSplitting, Vasilek.Collisions,
               Vasilek.VlasovPoisson1D1V, Vasilek.PoissonFourier1D, Vasilek.FDTD1D],
    checkdocs = :exports,
    doctest = false,
    repo = Documenter.Remotes.GitHub("korzhimanov", "Vasilek.jl"),
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true",
                             edit_link = "master"),
    pages = [
        "index.md",
        "Manual" => ["normalization.md", "driver-notes.md", "migration-0.2.md"],
        "Verification" => ["verification.md", ["studies/$stem.md" for stem in STUDIES]...],
        "api.md",
    ],
)

if EXECUTE
    deploydocs(; repo = "github.com/korzhimanov/Vasilek.jl.git", devbranch = "master")
else
    @info "The studies were not run, so the site is not deployed. Set VASILEK_DOCS_EXECUTE=1 to build the published form."
end
