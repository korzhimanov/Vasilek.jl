"""
    Collisions

Collision operators, as types, matching the convention of the advection
schemes in `Vasilek.Advection`: the operator is an immutable value and
`collide!` writes into an explicit destination.
"""
module Collisions

# One `workspace` for the package: a function of its own here was a second
# generic that the exported name never reached, so `workspace(BGK(τ), n)`
# was a MethodError.
import ..workspace

# `Landau1P` is defined here but not exported: it is a one-dimensional model of
# the Landau operator, with a bath at a temperature of its own; see its
# docstring.
export AbstractCollisionOperator, collide!, BGK

"""
    AbstractCollisionOperator

A velocity-space collision operator. Advance one step with

    collide!(dest, src, op, v, Δt[, ws])
"""
abstract type AbstractCollisionOperator end

"""
    workspace(op, n[, T])

Scratch for `op` at `n` velocity points, of the operator's element type unless
`T` is given, or `nothing`. `collide!` refuses one narrower than the data, so
Float64 data under `BGK(0.5f0)` needs `workspace(op, n, Float64)`; the
five-argument `collide!` builds its own in the data's type.
"""
workspace(::AbstractCollisionOperator, ::Integer, ::Type = Float64) = nothing

# The scratch follows the data, as the four-argument `advect!`'s does: one built
# from the operator's type would compute a Float64 line's Maxwellian in Float32
# under `BGK(0.5f0)`.
"""
    collide!(dest, src, op, v, Δt[, ws])

Advance the velocity line `src`, sampled at the nodes `v`, by a time `Δt` under
the collision operator `op`, and write the result into `dest`, which is
returned. `dest`, `src` and `v` have the same length.

`ws` is scratch from [`workspace`](@ref)`(op, length(v), float(eltype(src)))`,
one per task: its buffers of the lengths that call gives them, its type no
narrower than the data's, and sharing no memory with `dest` or `src`, which
[`_validate_collide`](@ref) checks. Without it a workspace of that type is
allocated for the call; pass one in any loop.
"""
collide!(dest, src, op::AbstractCollisionOperator, v, Δt) =
    collide!(dest, src, op, v, Δt, workspace(op, length(dest), float(eltype(src))))

"""
    _validate_collide(dest, src, op, v, ws)

The argument check at the top of every six-argument `collide!`. Each of these
gave a wrong answer, or an error from the wrong place, before it was checked:

  * `dest`, `src` and `v` of different lengths, or not one-based;
  * a workspace whose buffers are not the lengths `workspace(op, length(v))`
    gives them: a shorter one threw a `BoundsError` once the step was computed,
    after `Landau1P`'s O(n²) loop, and a longer `Landau1P` one ran;
  * a workspace narrower than the data, which the operator computes into: a
    Float32 one, as `workspace(BGK(0.5f0), n)` is, rounded a Float64 line's
    Maxwellian to Float32, 9.4e-9 off;
  * `dest` or `src` sharing memory with the workspace, which the operator
    writes while it reads `src` and before it writes `dest`;
  * `dest` sharing memory with `src` other than as `src` itself: `Landau1P`
    writes `dest[i]` before it reads `src[i+1]`, and `dest` one node on from
    `src` in the same array came back 1.005 off. `dest` may be `src` itself,
    and the step is taken in place.

An operator's workspace type has a method of `_scratch_lengths`,
`_scratch_eltype`, `_scratch_aliases` and `_err_scratch_alias`. The checks are
comparisons outside the loop and allocate nothing; the error paths are
`@noinline`.
"""
@inline function _validate_collide(dest, src, op::AbstractCollisionOperator, v, ws)
    Base.require_one_based_indexing(dest, src, v)
    n = length(src)
    length(dest) == n == length(v) || _err_collide_length(length(dest), n, length(v))
    _scratch_lengths(ws, n)
    T, S = _scratch_eltype(ws), float(eltype(src))
    promote_type(T, S) === T || _err_collide_narrow(T, S)
    _scratch_aliases(ws, dest, src) && _err_scratch_alias(ws)
    dest !== src && Base.mightalias(dest, src) && _err_collide_overlap()
    return nothing
end

_shares(a, dest, src) = Base.mightalias(a, dest) || Base.mightalias(a, src)

@noinline _err_collide_length(nd, ns, nv) = throw(DimensionMismatch(
    "collide! needs dest, src and v of one length, got $nd, $ns and $nv"))
@noinline _err_scratch_length(name, got, want) = throw(DimensionMismatch(
    "the workspace's $name has length $got, but this call needs $want; build it with " *
    "workspace(op, length(v), float(eltype(src)))"))
@noinline _err_collide_overlap() = throw(ArgumentError(
    "collide! requires dest to be src itself or not to share memory with it: the " *
    "operator writes dest while it still reads src's neighbours"))
@noinline _err_collide_narrow(T, S) = throw(ArgumentError(
    "the workspace is $T, narrower than the $S data: the operator computes into the " *
    "workspace's type and would round the data to it; build it with " *
    "workspace(op, length(v), $S)"))

include(joinpath("operators", "bgk.jl"))
include(joinpath("operators", "landau1p.jl"))

end # module
