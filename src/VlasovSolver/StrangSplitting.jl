module StrangSplitting

using ..Advection: AbstractAdvection1D, advect!
using LinearAlgebra: transpose!
import ..workspace

export strang_step!, make_time_step_2d!

"""
    StrangWorkspace

Scratch for [`strang_step!`](@ref): the transposed copy of `f` the second
direction sweeps over, a line buffer per direction, and each scheme's own
workspace. One per task.
"""
struct StrangWorkspace{T, WX, WV}
    ft::Matrix{T}
    bx::Vector{T}
    bv::Vector{T}
    wsx::WX
    wsv::WV
end

"""
    workspace(scheme_x, scheme_v, f::AbstractMatrix)

A [`StrangWorkspace`](@ref) for stepping `f` with the two schemes.
"""
function workspace(sx::AbstractAdvection1D, sv::AbstractAdvection1D, f::AbstractMatrix)
    T = float(eltype(f))
    nx, nv = size(f)
    return StrangWorkspace(Matrix{T}(undef, nv, nx), Vector{T}(undef, nx), Vector{T}(undef, nv),
                           workspace(sx, nx, T), workspace(sv, nv, T))
end

"""
    strang_step!(f, scheme_x, scheme_v, cx, cv, ws)

One Strang-split step `X(Δt/2) V(Δt) X(Δt/2)` of the matrix `f`, whose columns
are lines in the first direction (`f[:, j]` is the line at the `j`-th node of
the second one) and whose rows are lines in the second.

  * `cx` is the full-step argument of `advect!` for each column -- a Courant
    number for most schemes, a displacement for `PFCNonUniform` -- so
    `length(cx) == size(f, 2)`. Each half step takes `cx[j]/2`.
  * `cv(f)` returns the argument for each row, `size(f, 1)` of them, from `f`
    as it stands after the first half step: that is where the field is solved.

On return `f` holds the whole step; `ws.ft` is scratch and holds nothing
meaningful. `f` is updated in place.
"""
function strang_step!(f::AbstractMatrix, sx::AbstractAdvection1D, sv::AbstractAdvection1D,
                      cx::AbstractVector, cv, ws::StrangWorkspace)
    nx, nv = size(f)
    length(cx) == nv || throw(DimensionMismatch(
        "cx has $(length(cx)) entries, f has $nv columns"))
    _sweep!(f, sx, cx, 1//2, ws.bx, ws.wsx)
    c = cv(f)
    length(c) == nx || throw(DimensionMismatch(
        "cv(f) returned $(length(c)) entries, f has $nx rows"))
    transpose!(ws.ft, f)
    _sweep!(ws.ft, sv, c, 1, ws.bv, ws.wsv)
    transpose!(f, ws.ft)
    _sweep!(f, sx, cx, 1//2, ws.bx, ws.wsx)
    return f
end

function _sweep!(g, scheme, c, fraction, buf, ws)
    for j in axes(g, 2)
        line = view(g, :, j)
        advect!(buf, line, scheme, c[j]*fraction, ws)
        copyto!(line, buf)
    end
    return g
end

"""
    make_time_step_2d!(f, vΔt, advect!)

The 0.1 form: `f = (f₁, f₂)` holds the state twice, `f₂` the transpose of `f₁`;
`advect![d](line, α)` advances one line in place and `vΔt[d](other)` returns
the per-line arguments. Prefer [`strang_step!`](@ref), which takes schemes and
workspaces.

On return both arrays hold the whole step. `f₂` used to be left one half step
behind `f₁`. When `f₂` shares memory with `f₁` -- `f₂ = f₁'`, as the
verification harness passes -- the transposes are skipped, since there is
nothing to copy.
"""
function make_time_step_2d!(f, vΔt, advect!)
    shared = Base.mightalias(f[1], f[2])
    for (fi, α) in zip(eachcol(f[1]), vΔt[1](f[2]))
        advect![1](fi, α/2)
    end
    shared || transpose!(f[2], f[1])
    for (fi, α) in zip(eachcol(f[2]), vΔt[2](f[1]))
        advect![2](fi, α)
    end
    shared || transpose!(f[1], f[2])
    for (fi, α) in zip(eachcol(f[1]), vΔt[1](f[2]))
        advect![1](fi, α/2)
    end
    shared || transpose!(f[2], f[1])
    return f
end

end
