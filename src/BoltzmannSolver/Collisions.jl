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

# `Landau1P` is defined here but not exported: it is experimental, see its
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
`T` is given, or `nothing`.
"""
workspace(::AbstractCollisionOperator, ::Integer, ::Type = Float64) = nothing

"""
    collide!(dest, src, op, v, Δt[, ws])

Advance the velocity line `src`, sampled at the nodes `v`, by a time `Δt` under
the collision operator `op`, and write the result into `dest`, which is
returned. `dest`, `src` and `v` have the same length.

`ws` is scratch from [`workspace`](@ref)`(op, length(v), float(eltype(src)))`,
one per task, and must not share memory with `src`. Without it a workspace of
that type is allocated for the call; pass one in any loop.
"""
# The scratch follows the data, as the four-argument `advect!`'s does: one built
# from the operator's type would compute a Float64 line's Maxwellian in Float32
# under `BGK(0.5f0)`.
collide!(dest, src, op::AbstractCollisionOperator, v, Δt) =
    collide!(dest, src, op, v, Δt, workspace(op, length(dest), float(eltype(src))))

include(joinpath("operators", "bgk.jl"))
include(joinpath("operators", "landau1p.jl"))

end # module
