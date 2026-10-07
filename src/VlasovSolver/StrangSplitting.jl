module StrangSplitting

using ..Advection: AbstractAdvection1D, advect!
using ..Collisions: AbstractCollisionOperator, collide!
using LinearAlgebra: transpose!
import ..workspace

export strang_step!, Collide

"""
    StrangWorkspace

Scratch for [`strang_step!`](@ref): the transposed copy of `f` the second
direction sweeps over, a line buffer per direction, each scheme's own workspace
and the collision operator's. One per task.

`O` is the type of the operator it was built for, `Nothing` for none, and
`strang_step!` refuses a hook whose operator is of another kind.
"""
struct StrangWorkspace{S, T, WX, WV, WC, O}
    ft::Matrix{S}
    bx::Vector{T}
    bv::Vector{T}
    wsx::WX
    wsv::WV
    wsc::WC
end

"""
    workspace(scheme_x, scheme_v, f::AbstractMatrix, collisions = nothing, T = float(eltype(f)))

A [`StrangWorkspace`](@ref) for stepping `f` with the two schemes, and with
`collisions`, a collision operator, through a [`Collide`](@ref) hook.

The transposed copy holds `f`'s own element type, so the transposes copy it
exactly. The lines are worked in `T`: the buffers and the schemes' and the
operator's scratch take it, and each step rounds into `f` once. A wider `T`
than `f`'s -- `Float64` under Float32 `f`, as `vlasov_poisson` passes -- keeps a
flux sum from rounding a line past a scheme's bounds.
"""
function workspace(sx::AbstractAdvection1D, sv::AbstractAdvection1D, f::AbstractMatrix,
                   collisions::Union{Nothing, AbstractCollisionOperator} = nothing,
                   ::Type{T} = float(eltype(f))) where {T}
    nx, nv = size(f)
    wsx, wsv = workspace(sx, nx, T), workspace(sv, nv, T)
    wsc = collisions === nothing ? nothing : workspace(collisions, nv, T)
    return StrangWorkspace{eltype(f), T, typeof(wsx), typeof(wsv), typeof(wsc), typeof(collisions)}(
        Matrix{eltype(f)}(undef, nv, nx), Vector{T}(undef, nx), Vector{T}(undef, nv), wsx, wsv, wsc)
end

# The operator type a workspace was built for, and whether a hook's operator is
# of that kind: the same type, its parameters aside, as `BGK(0.5f0)`'s scratch
# serves `BGK(0.5)`. Decided from the types, so it costs nothing per step.
_built_for(::StrangWorkspace{S, T, WX, WV, WC, O}) where {S, T, WX, WV, WC, O} = O
_fits(ws::StrangWorkspace, op) = Base.typename(_built_for(ws)) === Base.typename(typeof(op))

@noinline _err_hook(ws, op) = throw(ArgumentError(
    "the hook collides with a $(nameof(typeof(op))), but the workspace was built for " *
    (_built_for(ws) === Nothing ? "no collision operator" : "a $(nameof(_built_for(ws)))") *
    "; build it with workspace(scheme_x, scheme_v, f, op)"))

"""
    Collide(op, v, Δt)

The collision hook of [`strang_step!`](@ref): every line of the second
direction, sampled at the nodes `v`, is collided under `op` for `Δt/2` on either
side of its advection, `C(Δt/2) V(Δt) C(Δt/2)`, so the split step stays
symmetric and second order. A value, built for each step; the operator's
scratch is the workspace's, from `workspace(scheme_x, scheme_v, f, op)`.
"""
struct Collide{O<:AbstractCollisionOperator, V<:AbstractVector, T}
    op::O
    v::V
    Δt::T
end

"""
    strang_step!(f, scheme_x, scheme_v, cx, cv, ws[, hook])

One Strang-split step `X(Δt/2) V(Δt) X(Δt/2)` of the matrix `f`, whose columns
are lines in the first direction (`f[:, j]` is the line at the `j`-th node of
the second one) and whose rows are lines in the second.

  * `cx` is the full-step argument of `advect!` for each column -- for an
    [`OnGrid`](@ref Vasilek.Advection.OnGrid) scheme a displacement, for a bare
    one its own fourth argument -- so `length(cx) == size(f, 2)`. Each half
    step takes `cx[j]/2`.
  * `cv(ft)` returns the argument for each line of the second direction,
    `size(f, 1)` of them. It is handed **`ft`, the state transposed**: `ft[:, i]`
    is the line at the `i`-th node of the first direction, as it stands after
    the first half step, which is where the field is solved.
  * `hook`, a [`Collide`](@ref), puts a collision operator on either side of
    the second sweep: `X(Δt/2) · C(Δt/2) V(Δt) C(Δt/2) · X(Δt/2)`. `ws` has to
    be built for that kind of operator, `workspace(scheme_x, scheme_v, f, op)`,
    or the call is an `ArgumentError` before `f` is touched.

The first direction is swept over the columns of `f` and the second over those
of `ws.ft`, so every line is contiguous. On return `f` holds the whole step;
`ws.ft` is scratch, and holds the state before the last half step. `f` is
updated in place.
"""
function strang_step!(f::AbstractMatrix, sx::AbstractAdvection1D, sv::AbstractAdvection1D,
                      cx::AbstractVector, cv, ws::StrangWorkspace,
                      hook::Union{Nothing, Collide} = nothing)
    nx, nv = size(f)
    length(cx) == nv || throw(DimensionMismatch(
        "cx has $(length(cx)) entries, f has $nv columns"))
    size(ws.ft) == (nv, nx) || throw(DimensionMismatch(
        "the workspace is for a $(reverse(size(ws.ft))) f, this one is $((nx, nv))"))
    if hook !== nothing
        length(hook.v) == nv || throw(DimensionMismatch(
            "the hook's v has $(length(hook.v)) nodes, f has $nv columns"))
        _fits(ws, hook.op) || _err_hook(ws, hook.op)
    end
    _sweep!(f, sx, cx, 1//2, ws.bx, ws.wsx)
    transpose!(ws.ft, f)
    c = cv(ws.ft)
    length(c) == nx || throw(DimensionMismatch(
        "cv(ft) returned $(length(c)) entries, f has $nx rows"))
    _kick!(ws.ft, sv, c, ws.bv, ws.wsv, hook, ws.wsc)
    transpose!(f, ws.ft)
    _sweep!(f, sx, cx, 1//2, ws.bx, ws.wsx)
    return f
end

# `c[j]*(1//2)` is `c[j]*0.5`, which is `c[j]/2` to the bit.
function _sweep!(g, scheme, c, fraction, buf, ws)
    for j in axes(g, 2)
        line = view(g, :, j)
        advect!(buf, line, scheme, c[j]*fraction, ws)
        copyto!(line, buf)
    end
    return g
end

function _kick!(g, scheme, c, buf, ws, ::Nothing, wsc)
    for j in axes(g, 2)
        line = view(g, :, j)
        advect!(buf, line, scheme, c[j], ws)
        copyto!(line, buf)
    end
    return g
end

function _kick!(g, scheme, c, buf, ws, hook::Collide, wsc)
    op, v, half = hook.op, hook.v, hook.Δt/2
    for j in axes(g, 2)
        line = view(g, :, j)
        collide!(buf, line, op, v, half, wsc)
        copyto!(line, buf)
        advect!(buf, line, scheme, c[j], ws)
        copyto!(line, buf)
        collide!(buf, line, op, v, half, wsc)
        copyto!(line, buf)
    end
    return g
end

end
