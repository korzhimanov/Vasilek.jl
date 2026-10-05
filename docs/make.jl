# Builds the documentation site: the README and the notes in docs/, the API
# reference, and the nine verification studies, rendered by Literate.jl with
# their output and figures. From the repository root:
#
#     julia --project=docs -e 'using Pkg; Pkg.instantiate()'     # once, Julia ≥ 1.11
#     julia --project=docs docs/make.jl                          # studies shown, not run
#     VASILEK_DOCS_EXECUTE=1 julia --project=docs docs/make.jl   # the published form
#
# and open docs/build/index.html.
#
# Running the studies is what puts their numbers and figures on the pages, and
# it costs the quarter of an hour Scripts.yml spends on them, so it is opt-in,
# as the extended test suite is through VASILEK_EXTENDED. Without it each study
# page shows its code under a note saying so, and the build still checks what a
# change can break: every docstring on the API page, every cross-reference
# resolving, every page parsing. `deploydocs` runs only after an executed build,
# so the published site never lacks its figures.

ENV["GKSwstype"] = "100"     # GR without a display, set before Plots loads it

using Documenter, Literate
using Vasilek

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC = joinpath(@__DIR__, "src")
const EXECUTE = get(ENV, "VASILEK_DOCS_EXECUTE", "0") == "1"
const BLOB = "https://github.com/korzhimanov/Vasilek.jl/blob/master/"

# Title on the site => script in verification/: the linear electrostatic cases,
# the instabilities, the kinetic and nonlinear ones, collisions, the
# electromagnetic one, and the advisory comparison.
const STUDIES = [
    "Landau damping"             => "landau-damping-1d1v",
    "Plasma oscillations"        => "plasma-oscillations-1d1v",
    "Two-stream instability"     => "two-stream",
    "Bump-on-tail instability"   => "bump-on-tail",
    "Plasma echo"                => "plasma-echo",
    "BGK equilibrium"            => "bgk-equilibrium",
    "Collisional Landau damping" => "collisional-damping",
    "Laser wakefield"            => "wakefield",
    "Scheme comparison"          => "scheme-comparison",
]

# ---------------------------------------------------------------- copied pages
# The README, the CHANGELOG and the notes in docs/ are linked by their paths
# from each other, from docstrings and from tests, so they stay where they are
# and are copied under src/ at build time. The copies are not committed.

# Where a file of the repository is on the site, relative to src/.
const SITE = Dict(
    "README.md"             => "index.md",
    "CHANGELOG.md"          => "changelog.md",
    "docs/verification.md"  => "verification.md",
    "docs/normalization.md" => "manual/normalization.md",
    "docs/migration-0.2.md" => "manual/migration-0.2.md",
    "docs/driver-notes.md"  => "manual/driver-notes.md",
    ("verification/$stem.jl" => "studies/$stem.md" for (_, stem) in STUDIES)...,
)

# The links of a page copied from `from` to `to`, rewritten for the site: a link
# to a file of the repository goes to the page that file became, or failing that
# to the file on GitHub, so that the same link works on both. Anything that does
# not name a file -- a URL, an anchor, an @ref, code that happens to read
# `](x)` -- is left as it is.
function relink(text, from, to)
    return replace(text, r"\]\(([^()\s]+)\)" => function (link)
        path, anchor = match(r"^([^#]*)(.*)$", chop(link; head = 2, tail = 1)).captures
        file = normpath(joinpath(dirname(from), path))
        (isempty(path) || occursin(r"^\w+:", path) || !ispath(file)) && return link
        repo = replace(relpath(file, ROOT), '\\' => '/')
        page = get(SITE, repo, nothing)
        new = page === nothing ? BLOB * repo :
              replace(relpath(joinpath(SRC, page), dirname(joinpath(SRC, to))), '\\' => '/')
        return "](" * new * anchor * ")"
    end)
end

# `EditURL` points the page's "Edit on GitHub" link at the original. It is
# relative to the copy, and Documenter refuses one that names no file.
function copy_page(from, to; meta = "")
    dest = joinpath(SRC, to)
    mkpath(dirname(dest))
    editurl = replace(relpath(from, dirname(dest)), '\\' => '/')
    write(dest, "```@meta\nEditURL = \"$editurl\"\n$meta```\n\n",
          relink(read(from, String), from, to))
    return to
end

copy_page(joinpath(ROOT, "README.md"), "index.md")
copy_page(joinpath(ROOT, "CHANGELOG.md"), "changelog.md")
copy_page(joinpath(@__DIR__, "verification.md"), "verification.md")
for name in ("normalization", "migration-0.2", "driver-notes")
    # The driver notes link its internals with @ref, which resolve in its module.
    meta = name == "driver-notes" ? "CurrentModule = Vasilek.VlasovPoisson1D1V\n" : ""
    copy_page(joinpath(@__DIR__, "$name.md"), "manual/$name.md"; meta)
end

# -------------------------------------------------------------------- studies
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

for (_, stem) in STUDIES
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
    # Every module with docstrings: `@docs workspace` gathers the methods'
    # docstrings from these alone, and `checkdocs` checks these alone.
    modules = [Vasilek, Vasilek.Advection, Vasilek.StrangSplitting,
               Vasilek.VlasovPoisson1D1V, Vasilek.PoissonFourier1D,
               Vasilek.FDTD1D, Vasilek.Collisions],
    repo = Documenter.Remotes.GitHub("korzhimanov", "Vasilek.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "master",
        # 150 kB of Markdown, over the HTML size limit; every other page is small.
        size_threshold_ignore = ["changelog.md"],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Normalization" => "manual/normalization.md",
            "The 1D1V driver" => "manual/driver-notes.md",
            "Migrating to 0.2" => "manual/migration-0.2.md",
        ],
        "Verification" => ["Overview" => "verification.md",
                           [title => "studies/$stem.md" for (title, stem) in STUDIES]...],
        "API reference" => "api.md",
        "Changelog" => "changelog.md",
    ],
    checkdocs = :all,
    warnonly = false,
)

if EXECUTE
    deploydocs(; repo = "github.com/korzhimanov/Vasilek.jl.git", devbranch = "master")
else
    @info "The studies were not run, so the site is not deployed. Set VASILEK_DOCS_EXECUTE=1 to build the published form."
end
