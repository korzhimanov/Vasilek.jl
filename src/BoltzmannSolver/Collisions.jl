"""
    Collisions

Collision operators, as types, matching the convention used by
[`Advection`](@ref): the operator is an immutable value and `collide!` writes
into an explicit destination.
"""
module Collisions

using NumericalIntegration

export AbstractCollisionOperator, collide!, BGK, Landau1P

"""
    AbstractCollisionOperator

A velocity-space collision operator. Advance one step with

    collide!(dest, src, op, v, Δt[, ws])
"""
abstract type AbstractCollisionOperator end

"""
    workspace(op, n)

Scratch for `op` at `n` velocity points, or `nothing`.
"""
workspace(::AbstractCollisionOperator, ::Integer) = nothing

collide!(dest, src, op::AbstractCollisionOperator, v, Δt) =
    collide!(dest, src, op, v, Δt, workspace(op, length(dest)))

include(joinpath("operators", "bgk.jl"))
include(joinpath("operators", "landau1p.jl"))

end # module
