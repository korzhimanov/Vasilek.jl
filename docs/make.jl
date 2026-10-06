# Builds the documentation into docs/build. Not deployed anywhere yet; CI builds
# it so that a broken reference or a docstring left out fails a pull request.
#
#     julia --project=docs -e 'using Pkg; Pkg.instantiate()'   # once, Julia ≥ 1.11
#     julia --project=docs docs/make.jl
using Documenter
using Vasilek

makedocs(;
    sitename = "Vasilek.jl",
    modules = [Vasilek, Vasilek.Advection, Vasilek.StrangSplitting, Vasilek.Collisions,
               Vasilek.VlasovPoisson1D1V, Vasilek.PoissonFourier1D, Vasilek.FDTD1D],
    checkdocs = :exports,
    doctest = false,
    pages = ["index.md", "normalization.md", "driver-notes.md", "migration-0.2.md", "api.md"],
)
