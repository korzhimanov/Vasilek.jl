"""
    BGK(τ; conservative = true)

Relaxation towards the local Maxwellian on timescale `τ`: the exact solution of
`∂f/∂t = (M − f)/τ` with `M` held fixed over the step,

    dest = src·exp(−Δt/τ) + (1 − exp(−Δt/τ))·M.

**`conservative = true` (the default)** takes `M` as the *discrete* Maxwellian
of the line (Mieussens, Math. Models Methods Appl. Sci. 10, 1121 (2000)): the
function `exp(a + bv + cv²)` whose density, momentum and energy, summed on the
grid with the cell widths the advection schemes conserve,

    Σᵢ wᵢ fᵢ (1, vᵢ, vᵢ²),   w = (v₂ − v₁, (v₃ − v₁)/2, …, vₙ − vₙ₋₁),

equal those of `src`. It is found by Newton's method on the convex dual,
starting from the continuous Maxwellian, so `collide!` conserves all three to
round-off on any grid and any window, and a sampled Maxwellian is a fixed point.

**`conservative = false`** samples the continuous Maxwellian of the trapezoid
moments, as the operator did before 0.2. It conserves only as far as the window
holds the Maxwellian's tails: a line relaxed on ±4 lost 0.8% of its density
and 7% of its energy where ±8 lost 1e-6.

An empty line (`n ≤ 0`) and one colder than the grid resolves (`T < Δv²`,
which includes a single-node spike) are returned unchanged. So is a line whose
moments no `exp(a + bv + cv²)` on the grid can match, where Newton does not
converge: all of the mass on the two end nodes is one, with `T` at the largest
the window allows. Leaving one such line alone conserves everything; throwing
would end a whole `vlasov_poisson` run over a single column.
"""
struct BGK{T<:AbstractFloat} <: AbstractCollisionOperator
    τ::T
    conservative::Bool
end
BGK(τ; conservative::Bool = true) = BGK(float(τ), conservative)

struct BGKWorkspace{T}
    maxwellian::Vector{T}
    moment::Vector{T}
end
workspace(::BGK{T}, n::Integer, ::Type{S} = T) where {T, S} =
    BGKWorkspace(Vector{S}(undef, n), Vector{S}(undef, n))

function collide!(dest, src, op::BGK, v, Δt, ws::BGKWorkspace)
    length(dest) == length(src) == length(v) || throw(DimensionMismatch(
        "collide! needs dest, src and v of one length, got $(length(dest)), " *
        "$(length(src)) and $(length(v))"))
    # The Maxwellian, and the sampled operator's moments, are written before
    # `src` is read for the last time: a `src` sharing their memory came back as
    # the Maxwellian alone, 0.08 from the step, or 0.64 off through `moment`.
    (Base.mightalias(src, ws.maxwellian) || Base.mightalias(src, ws.moment)) &&
        _err_scratch_alias()
    e = exp(-Δt/op.τ)
    M = ws.maxwellian
    ok = op.conservative ? _discrete_maxwellian!(M, src, v) :
                           _sampled_maxwellian!(M, ws.moment, src, v)
    ok || return copyto!(dest, src)
    @. dest = src*e + (1 - e)*M
    return dest
end

@noinline _err_scratch_alias() = throw(ArgumentError(
    "collide! requires src not to share memory with the workspace: the operator " *
    "writes its Maxwellian and moments there before it is done reading src"))

# The coarsest spacing: a line whose temperature is below its square is not
# resolved, and its Maxwellian would be a spike the grid cannot integrate.
_coarsest(v) = maximum(v[i+1] - v[i] for i in 1:length(v)-1)

# Cell width of node i, the weight of the flux-form sums.
@inline function _width(v, i)
    n = length(v)
    return i == 1 ? v[2] - v[1] : i == n ? v[n] - v[n-1] : (v[i+1] - v[i-1])/2
end

function _trapezoid(v, y)
    s = zero(promote_type(eltype(v), eltype(y)))
    for i in 1:length(v)-1
        s += (v[i+1] - v[i])*(y[i] + y[i+1])
    end
    return s/2
end

"The 0.1 operator's Maxwellian: continuous, from trapezoid moments."
function _sampled_maxwellian!(M, moment, src, v)
    n = _trapezoid(v, src)
    n > 0 || return false
    @. moment = v*src
    u = _trapezoid(v, moment)/n
    @. moment = (v - u)^2 * src
    T = _trapezoid(v, moment)/n
    T ≥ _coarsest(v)^2 || return false
    @. M = n/sqrt(2π*T)*exp(-(v - u)^2/(2T))
    return true
end

"""
    _discrete_maxwellian!(M, f, v)

Write into `M` the discrete Maxwellian of `f` on `v`, or return `false` when the
line has none to relax to. In the variable `ξ = (v − u)/√T` the target is
`M = exp(α₀ + α₁ξ + α₂ξ²)` with `Σ w M (1, ξ, ξ²) = Σ w f (1, ξ, ξ²) = ρ`, the
minimum of the convex `Φ(α) = Σ w exp(α·φ) − α·ρ`. Newton with a backtracking
line search, from the continuous Maxwellian `α = (log(n/√(2πT)), 0, −1/2)`;
on a resolved line it takes two or three steps.
"""
function _discrete_maxwellian!(M, f, v)
    R = float(promote_type(eltype(f), eltype(v)))
    n = zero(R); p = zero(R)
    for i in eachindex(v)
        w = _width(v, i)
        n += w*f[i]; p += w*f[i]*v[i]
    end
    n > 0 || return false
    u = p/n
    q = zero(R)
    for i in eachindex(v)
        q += _width(v, i)*f[i]*(v[i] - u)^2
    end
    T = q/n
    T ≥ _coarsest(v)^2 || return false
    s = sqrt(T)

    ρ₀ = n; ρ₁ = zero(R); ρ₂ = zero(R)
    for i in eachindex(v)
        ξ = (v[i] - u)/s
        wf = _width(v, i)*f[i]
        ρ₁ += wf*ξ; ρ₂ += wf*ξ^2
    end

    α₀, α₁, α₂ = log(n/sqrt(2π*T)), zero(R), -one(R)/2
    tol = 64*eps(R)*n
    for _ in 1:50
        m = _maxwellian_moments(v, u, s, α₀, α₁, α₂)
        g₀, g₁, g₂ = m[1] - ρ₀, m[2] - ρ₁, m[3] - ρ₂
        max(abs(g₀), abs(g₁), abs(g₂)) ≤ tol && break
        # Newton direction: H d = g, H the symmetric moment matrix (m₀ … m₄)
        d₀, d₁, d₂ = _solve3(m[1], m[2], m[3], m[4], m[5], g₀, g₁, g₂)
        Φ = m[1] - (α₀*ρ₀ + α₁*ρ₁ + α₂*ρ₂)
        t = one(R)
        while true
            β₀, β₁, β₂ = α₀ - t*d₀, α₁ - t*d₁, α₂ - t*d₂
            Φt = _maxwellian_moments(v, u, s, β₀, β₁, β₂)[1] - (β₀*ρ₀ + β₁*ρ₁ + β₂*ρ₂)
            # Near the optimum Φ changes by less than its own round-off, so a
            # strict decrease would refuse the full Newton step at random and
            # stall; within that round-off the step is taken.
            if Φt ≤ Φ + 64*eps(R)*(abs(Φ) + n) || t < 1e-6
                α₀, α₁, α₂ = β₀, β₁, β₂
                break
            end
            t /= 2
        end
    end
    m = _maxwellian_moments(v, u, s, α₀, α₁, α₂)
    # Not realisable on this grid: no Maxwellian to relax to, as for n ≤ 0.
    max(abs(m[1] - ρ₀), abs(m[2] - ρ₁), abs(m[3] - ρ₂)) ≤ 1024*eps(R)*n || return false
    for i in eachindex(v)
        ξ = (v[i] - u)/s
        M[i] = exp(α₀ + ξ*(α₁ + ξ*α₂))
    end
    return true
end

"`Σ w M ξᵏ` for k = 0…4, M = exp(α₀ + α₁ξ + α₂ξ²)."
function _maxwellian_moments(v, u, s, α₀, α₁, α₂)
    m₀ = m₁ = m₂ = m₃ = m₄ = zero(α₀)
    for i in eachindex(v)
        ξ = (v[i] - u)/s
        wM = _width(v, i)*exp(α₀ + ξ*(α₁ + ξ*α₂))
        m₀ += wM; m₁ += wM*ξ; m₂ += wM*ξ^2; m₃ += wM*ξ^3; m₄ += wM*ξ^4
    end
    return (m₀, m₁, m₂, m₃, m₄)
end

# Solve [a b c; b c d; c d e] x = r by Cramer's rule; the matrix is a Hankel
# moment matrix, symmetric positive definite for a positive M.
function _solve3(a, b, c, d, e, r₀, r₁, r₂)
    det = a*(c*e - d*d) - b*(b*e - c*d) + c*(b*d - c*c)
    x₀ = (r₀*(c*e - d*d) - b*(r₁*e - d*r₂) + c*(r₁*d - c*r₂))/det
    x₁ = (a*(r₁*e - d*r₂) - r₀*(b*e - c*d) + c*(b*r₂ - r₁*c))/det
    x₂ = (a*(c*r₂ - r₁*d) - b*(b*r₂ - r₁*c) + r₀*(b*d - c*c))/det
    return x₀, x₁, x₂
end
