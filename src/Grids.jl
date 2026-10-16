"""
    Grids

One-dimensional grids, as types: a [`UniformGrid`](@ref), whose cells share one
spacing, and a [`NonUniformGrid`](@ref), whose cells each have a width of their
own. Whether a grid is uniform is decided once, where it is built, and the code
that needs a spacing dispatches on the type: [`OnGrid`](@ref
Vasilek.Advection.OnGrid) takes a scheme that wants a Courant number on a
`UniformGrid` only, and the driver's field solve takes its `x` grid the same way.

A range is a uniform grid as it stands, its `step` the spacing to the bit. A
vector of nodes is classified by [`grid`](@ref), the one place that compares
spacings, to the rounding the nodes carry.
"""
module Grids

export AbstractGrid1D, UniformGrid, NonUniformGrid, grid,
       nodes, hasnodes, spacing, widths, narrowest, cell_widths, samecells

"""
    AbstractGrid1D{T}

A one-dimensional grid whose spacings or widths are of the floating-point type
`T`: a [`UniformGrid`](@ref) or a [`NonUniformGrid`](@ref). `length(g)` is the
number of its cells, one per node.
"""
abstract type AbstractGrid1D{T<:AbstractFloat} end

"""
    UniformGrid(r::AbstractRange)
    UniformGrid(z::AbstractVector)
    UniformGrid(h, nodes)
    UniformGrid(h, n)

A grid whose `n` cells share the spacing `h`, and optionally the nodes it was
built from. `h` is *the* spacing: what a displacement is divided by for a
Courant number, the field solve's sampling interval, and every cell's width.
The nodes agree with it to their own rounding only.

  * From a range, `h` is `step(r)`, exactly the step it was written with:
    `UniformGrid(-4.0:0.2:4.0)` has `h = 0.2`, where the collected nodes' first
    difference is `0.20000000000000018`. Its nodes, collected in `T`, are where
    `f₀` is sampled, and are checked as a vector's are: a Float16 range over
    [−4, 4] whose nodes round by 1% of a cell is not uniform.
  * From a vector of nodes, the grid has to be uniform to rounding, as
    [`grid`](@ref) decides, or the call is an `ArgumentError` naming the first
    node that is not; `h` is then the mean spacing, `(z[end] - z[1])/(n - 1)`,
    so that `n` cells of it span the nodes' periodic box to rounding. The first
    difference carries the rounding of the largest node: `n` cells of it were
    2.3e-10 off the box of `1e6 .+ (0:63)*0.1`, the mean 7.4e-12.
  * `UniformGrid(h, nodes)` takes the spacing as given, for nodes uniform to
    rounding that it is the spacing of, to their rounding: cell centres
    `(j - 1/2)L/N` with `h = L/N`, whose differences round otherwise. An `h`
    the nodes do not have is an `ArgumentError`, since the field solve and the
    Courant numbers would take it while the sums take the nodes' cells.
  * `UniformGrid(h, n)` has no nodes: enough for a scheme on a grid, not for
    the driver, which samples `f₀` and the ions at the nodes.

`T` is `float` of the nodes' type and `h`'s: a Float32 grid stays Float32, a
Rational or an integer one becomes Float64.
"""
struct UniformGrid{T<:AbstractFloat, X<:Union{Nothing, Vector{T}}} <: AbstractGrid1D{T}
    h::T
    n::Int
    nodes::X
    # Every UniformGrid is built through here, and checked in `T`, as it holds
    # its spacing and nodes: a spacing can round to zero in `T`, nodes to equal
    # ones.
    function UniformGrid{T,X}(h, n, nodes) where {T<:AbstractFloat, X<:Union{Nothing, Vector{T}}}
        h = T(h)
        isfinite(h) && h > 0 || throw(ArgumentError(
            "a UniformGrid needs a finite, positive spacing, got h = $h"))
        n ≥ 1 || throw(ArgumentError("a UniformGrid needs at least one cell, got n = $n"))
        if nodes !== nothing
            nodes = convert(Vector{T}, nodes)
            length(nodes) == n || throw(DimensionMismatch(
                "a UniformGrid of $n cells given $(length(nodes)) nodes"))
            _check_nodes(nodes)
        end
        return new{T,X}(h, n, nodes)
    end
end

function UniformGrid(h::Real, n::Integer)
    T = float(typeof(h))
    return UniformGrid{T, Nothing}(T(h), n, nothing)
end

function UniformGrid(r::AbstractRange)
    x = _range_nodes(r)
    h, j = _uniform_spacing(x)
    j == 0 || _err_nonuniform(x, h, j)
    return _uniform_range(r, x)
end

# A range's nodes, collected in its float type, which is where they are used.
function _range_nodes(r::AbstractRange)
    _check_nodes(r)
    return collect(float(eltype(r)), r)
end
_uniform_range(r, x::Vector{T}) where {T} = UniformGrid{T, Vector{T}}(T(step(r)), length(r), x)

function UniformGrid(h::Real, z::AbstractVector)
    _check_nodes(z)
    hz, j = _uniform_spacing(z)
    j == 0 || _err_nonuniform(z, hz, j)
    T = float(promote_type(typeof(h), eltype(z)))
    m = _mean_spacing(z)
    _same_width(T(h), T(m), _rtol(T), maximum(abs, z)) || throw(ArgumentError(
        "h = $h is not the spacing of the nodes, $m on average, to their rounding"))
    return _uniform(h, z)
end

function UniformGrid(z::AbstractVector)
    _check_nodes(z)
    h, j = _uniform_spacing(z)
    j == 0 || _err_nonuniform(z, h, j)
    return _uniform(_mean_spacing(z), z)
end

# Nodes already checked uniform, and their spacing.
function _uniform(h, z)
    T = float(promote_type(typeof(h), eltype(z)))
    return UniformGrid{T, Vector{T}}(T(h), length(z), collect(T, z))
end

# In the nodes' own arithmetic: exact for Rational nodes.
_mean_spacing(z) = (z[end] - z[1])/(length(z) - 1)

"""
    NonUniformGrid(z::AbstractVector)
    NonUniformGrid(; widths)

A grid whose cells each have a width of their own: from nodes, the widths
[`cell_widths`](@ref) gives them, computed once; or from the widths alone, as a
[`PFCNonUniform`](@ref Vasilek.Advection.PFCNonUniform) carries them. `hmin`
is the narrowest.

Nothing is checked about uniformity: a `NonUniformGrid` of equal widths is
still non-uniform to every method that dispatches on it, and a scheme that
takes a Courant number is refused on it.
"""
struct NonUniformGrid{T<:AbstractFloat, X<:Union{Nothing, Vector{T}}} <: AbstractGrid1D{T}
    widths::Vector{T}
    hmin::T
    nodes::X
    # Checked in `T`, as `UniformGrid` is. `convert` keeps a `Vector{T}` as it
    # is: `OnGrid(p::PFCNonUniform)` shares the scheme's widths.
    function NonUniformGrid{T,X}(widths, hmin, nodes) where {T<:AbstractFloat, X<:Union{Nothing, Vector{T}}}
        widths, hmin = convert(Vector{T}, widths), T(hmin)
        isempty(widths) && throw(ArgumentError("a NonUniformGrid needs at least one cell"))
        all(w -> isfinite(w) && w > 0, widths) || throw(ArgumentError(
            "a NonUniformGrid needs finite, positive cell widths; got extrema $(extrema(widths))"))
        hmin == minimum(widths) || throw(ArgumentError(
            "NonUniformGrid: hmin = $hmin is not the narrowest width, $(minimum(widths))"))
        if nodes !== nothing
            nodes = convert(Vector{T}, nodes)
            length(nodes) == length(widths) || throw(DimensionMismatch(
                "a NonUniformGrid of $(length(widths)) cells given $(length(nodes)) nodes"))
            _check_nodes(nodes)
        end
        return new{T,X}(widths, hmin, nodes)
    end
end

function NonUniformGrid(z::AbstractVector)
    _check_nodes(z)
    T = float(eltype(z))
    x = collect(T, z)
    w = cell_widths(x)
    return NonUniformGrid{T, Vector{T}}(w, minimum(w), x)
end

function NonUniformGrid(; widths::AbstractVector{<:Real})
    T = float(eltype(widths))
    w = collect(T, widths)
    isempty(w) && throw(ArgumentError("a NonUniformGrid needs at least one cell"))
    return NonUniformGrid{T, Nothing}(w, minimum(w), nothing)
end

# --------------------------------------------------------------- the one rule

"""
    CAP

The most the ulp allowance of [`_within_rounding`](@ref) may be, as a fraction
of a cell: 0.5%. Spacings carry the rounding of the coordinates they are
differences of, which grows with the coordinates' magnitude rather than with the
spacing; past this fraction a grid is not uniform whatever its rounding, since
a Courant number taken from its first spacing would be off by as much
elsewhere. Without it 8 ulps of a Float16 grid's largest node were 16% of a
cell, and a Float16 grid whose cells differ by 12% ran `Upwind` at one Courant
number. The acceptance it gives is measured in `docs/src/driver-notes.md`,
"Grids".
"""
const CAP = 5e-3

"""
    _within_rounding(d, h, scale)

Whether the spacing `d` is `h` to rounding: `1e-10` of `h`, or 8 ulps of `scale`
in `scale`'s own type -- the largest coordinate the spacings are differences of
-- but never more than [`CAP`](@ref) of `h`. NaN is never within rounding.
"""
_within_rounding(d, h, scale) =
    d == h || abs(d - h) ≤ max(1e-10*abs(h), min(8eps(float(scale)), CAP*abs(h)))

"""
    _uniform_spacing(z) -> (h, j)

The first spacing of the nodes `z`, `z[2] - z[1]` in their own arithmetic, and
the first node `j` whose spacing is not that to rounding,
[`_within_rounding`](@ref) at the scale of the largest node; `j = 0` when every
one is.
"""
function _uniform_spacing(z::AbstractVector)
    h = z[2] - z[1]
    scale = maximum(abs, z)
    for j in 2:length(z)-1
        _within_rounding(z[j+1] - z[j], h, scale) || return (h, j)
    end
    return (h, 0)
end

"""
    grid(z)

The grid `z` is: itself for an [`AbstractGrid1D`](@ref), and for a vector of
nodes or a range a [`UniformGrid`](@ref) if every spacing is the first to
rounding and a [`NonUniformGrid`](@ref) otherwise. A vector's spacing is its
mean; a range's is its `step`, and its nodes are checked, collected in its
float type, since they are where the run samples its data. This is where the
package decides that a vector is uniform, and the only place.

Rounding is `1e-10` of the first spacing, or 8 ulps of the largest node in the
nodes' own type, whichever is more, and never more than 0.5% of a cell: the
spacings of a collected range scatter by a few ulps of its largest node, which
for Float32 nodes up to 25 is 2.8e-6 of a π/8 spacing and for Float64 ones near
1e6 is 1.2e-9 of 0.1, and both are uniform; a grid stretched by a millionth of a
cell, a Float16 grid whose spacings differ by 1%, and a Float32 grid that
reaches 10⁵ cells from the origin are not. The measurements are in `docs/src/driver-notes.md`, "Grids".

Nodes have to be finite and increase strictly, at least two of them. A vector
of cell widths is not nodes: widths that do not increase are refused as such,
but widths that do, as a graded grid's, are taken for nodes like any others.
"""
grid(g::AbstractGrid1D) = g
function grid(r::AbstractRange)
    x = _range_nodes(r)
    h, j = _uniform_spacing(x)
    return j == 0 ? _uniform_range(r, x) : NonUniformGrid(x)
end
function grid(z::AbstractVector)
    _check_nodes(z)
    h, j = _uniform_spacing(z)
    return j == 0 ? _uniform(_mean_spacing(z), z) : NonUniformGrid(z)
end

function _check_nodes(z)
    Base.require_one_based_indexing(z)
    n = length(z)
    n ≥ 2 || throw(ArgumentError("a grid needs at least 2 nodes, got $n"))
    for j in 1:n
        isfinite(z[j]) || throw(ArgumentError("a grid's nodes must be finite, node $j is $(z[j])"))
    end
    for j in 1:n-1
        z[j+1] > z[j] || _err_increasing(z, j)
    end
    return nothing
end

@noinline _err_increasing(z::AbstractRange, j) = throw(ArgumentError(
    "a range given for a grid has to increase; this one steps by $(step(z))"))
@noinline _err_increasing(z, j) = throw(ArgumentError(
    "a grid's nodes must increase, node $(j+1) is $(z[j+1]) after $(z[j]). " *
    "A vector given for a grid is its nodes, not its cell widths: give widths " *
    "as NonUniformGrid(; widths), or a uniform spacing as (h, n)"))

@noinline _err_nonuniform(z, h, j) = throw(ArgumentError(
    "the nodes are not uniformly spaced: spacing $(z[j+1] - z[j]) at node $j, " *
    "$h at node 1; NonUniformGrid takes them, or grid(z) decides"))

# ------------------------------------------------------------------ accessors

Base.length(g::UniformGrid) = g.n
Base.length(g::NonUniformGrid) = length(g.widths)
Base.eltype(::Type{<:AbstractGrid1D{T}}) where {T} = T

"""
    hasnodes(g)

Whether the grid `g` was built from its nodes, and holds them.
"""
hasnodes(g::AbstractGrid1D) = g.nodes !== nothing

"""
    nodes(g)

The nodes the grid `g` was built from, in its own type; an `ArgumentError` for
a grid built from its spacing or its widths alone.
"""
nodes(g::AbstractGrid1D) = hasnodes(g) ? g.nodes : _err_nodes(g)

@noinline _err_nodes(g) = throw(ArgumentError(
    "$g has no nodes: it was built from its spacing or its widths alone"))

"""
    spacing(g::UniformGrid)

The spacing of a uniform grid. A [`NonUniformGrid`](@ref) has none, and no
method.
"""
spacing(g::UniformGrid) = g.h

"""
    widths(g)

The width of each cell: `h` for every cell of a [`UniformGrid`](@ref), a vector
built on each call; the widths a [`NonUniformGrid`](@ref) holds.
"""
widths(g::UniformGrid) = fill(g.h, g.n)
widths(g::NonUniformGrid) = g.widths

"""
    narrowest(g)

The narrowest cell of the grid `g`: its spacing if uniform.
"""
narrowest(g::UniformGrid) = g.h
narrowest(g::NonUniformGrid) = g.hmin

"""
    cell_widths(z)
    cell_widths(g::AbstractGrid1D)

The width each node of the grid `z` stands for: the spacing at the two ends,
and half the distance between its neighbours inside. These are the weights of
the sums `Σ f ΔvΔx` that the flux-form schemes conserve. A floating-point grid
gets widths of its own element type. Of a grid, its [`widths`](@ref); of a
range, its step in every cell, as of its `UniformGrid` -- not the differences of
its nodes, which round in the last bit.
"""
cell_widths(z::AbstractVector) = vcat([z[2]-z[1]], (z[3:end] - z[1:end-2])/2, [z[end]-z[end-1]])
cell_widths(r::AbstractRange) = widths(UniformGrid(r))
cell_widths(g::AbstractGrid1D) = widths(g)

# ----------------------------------------------------------- the same grid

# Spacings of two types compare to 4 ulps of the coarser, as `Yee1D` compares
# its steps: an `OnGrid` built on a Float32 grid's spacing is that grid in
# Float64 too.
_rtol(::Type{T}) where {T<:AbstractFloat} = max(1e-12, 4eps(T))
_rtol(a::AbstractGrid1D, b::AbstractGrid1D) = max(_rtol(eltype(typeof(a))), _rtol(eltype(typeof(b))))

# The rounding a grid's spacing carries: that of its largest node, in the
# nodes' own type, for a grid built from nodes; none but `1e-10` of the spacing
# for one built from its spacing or widths.
_scale(g::AbstractGrid1D) = hasnodes(g) ? maximum(abs, g.nodes) : zero(eltype(typeof(g)))

# Two widths or spacings of the same grid: to the coarser type's 4 ulps, or to
# the rounding of the nodes of either, each in its own type -- a Float32 grid's
# taken in Float64 lost its rounding.
_same_width(d, w, rtol, scales...) =
    isapprox(d, w; rtol) || any(s -> _within_rounding(d, w, s), scales)

"""
    samecells(a, b)

Whether the grids `a` and `b` have the same cells: as many, and the same
spacing, or the same widths cell by cell, to the rounding either carries -- 4
ulps of the coarser type, or the rounding of the nodes they were built from,
by the rule [`grid`](@ref) applies. A uniform grid and a non-uniform one never
have.
"""
samecells(a::UniformGrid, b::UniformGrid) =
    a.n == b.n && _same_width(a.h, b.h, _rtol(a, b), _scale(a), _scale(b))
function samecells(a::NonUniformGrid, b::NonUniformGrid)
    length(a) == length(b) || return false
    rtol, sa, sb = _rtol(a, b), _scale(a), _scale(b)
    return all(_same_width(p, q, rtol, sa, sb) for (p, q) in zip(a.widths, b.widths))
end
samecells(::AbstractGrid1D, ::AbstractGrid1D) = false

function Base.show(io::IO, g::UniformGrid{T}) where {T}
    print(io, "UniformGrid{", T, "}(h = ", g.h, ", n = ", g.n)
    hasnodes(g) && print(io, ", nodes ", first(g.nodes), "…", last(g.nodes))
    print(io, ")")
end
function Base.show(io::IO, g::NonUniformGrid{T}) where {T}
    print(io, "NonUniformGrid{", T, "}(n = ", length(g), ", widths ", g.hmin, "…",
          maximum(g.widths))
    hasnodes(g) && print(io, ", nodes ", first(g.nodes), "…", last(g.nodes))
    print(io, ")")
end

end # module
