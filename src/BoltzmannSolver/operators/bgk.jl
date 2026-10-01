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
    # An empty line (vacuum) has no drift or temperature to relax towards, and
    # 0/0 would fill it with NaN; nor does one colder than the grid resolves.
    # Collisions leave both as they are. "Resolves" means T ≥ Δv²: a single-node
    # spike has T = 0 only up to round-off -- off a dyadic node u misses v by an
    # ulp and T comes out near 1e-34, whose Maxwellian is a 1e14 spike -- and
    # below Δv² the sampled Maxwellian no longer integrates to n (at T = Δv² the
    # trapezoid holds it to 5e-9).
    n > 0 || return copyto!(dest, src)
    @. moment = v*src
    u = integrate(v, moment)/n
    @. moment = (v - u)^2 * src
    T = integrate(v, moment)/n
    Δv = maximum(v[i+1] - v[i] for i in 1:length(v)-1)
    T ≥ Δv^2 || return copyto!(dest, src)

    @. M = n/sqrt(2π*T)*exp(-(v - u)^2/(2T))
    @. dest = src*e + (1.0 - e)*M
    return dest
end
