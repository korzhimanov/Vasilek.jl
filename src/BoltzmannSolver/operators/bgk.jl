"""
    BGK(τ)

Relaxation towards the local Maxwellian on timescale `τ`.
"""
struct BGK{T<:AbstractFloat} <: AbstractCollisionOperator
    τ::T
end
BGK(τ) = BGK(float(τ))

struct BGKWorkspace{T}
    maxwellian::Vector{T}
    moment::Vector{T}
end
workspace(::BGK{T}, n::Integer, ::Type{S} = T) where {T, S} =
    BGKWorkspace(Vector{S}(undef, n), Vector{S}(undef, n))

function collide!(dest, src, op::BGK, v, Δt, ws::BGKWorkspace)
    e = exp(-Δt/op.τ)
    M, moment = ws.maxwellian, ws.moment

    n = integrate(v, src)
    # An empty line (vacuum) has no drift or temperature to relax towards, and
    # 0/0 would fill it with NaN; nor does one colder than the grid resolves.
    # Collisions leave both as they are.
    n > 0 || return copyto!(dest, src)
    @. moment = v*src
    u = integrate(v, moment)/n
    @. moment = (v - u)^2 * src
    T = integrate(v, moment)/n
    T > 0 || return copyto!(dest, src)

    @. M = n/sqrt(2π*T)*exp(-(v - u)^2/(2T))
    @. dest = src*e + (1 - e)*M
    return dest
end
