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
workspace(::BGK, n::Integer) = BGKWorkspace(Vector{Float64}(undef, n), Vector{Float64}(undef, n))

function collide!(dest, src, op::BGK, v, Δt, ws::BGKWorkspace)
    e = exp(-Δt/op.τ)
    M, moment = ws.maxwellian, ws.moment

    n = integrate(v, src)
    @. moment = v*src
    u = integrate(v, moment)/n
    @. moment = (v - u)^2 * src
    T = integrate(v, moment)/n

    @. M = n/sqrt(2π*T)*exp(-(v - u)^2/(2T))
    @. dest = src*e + (1.0 - e)*M
    return dest
end
