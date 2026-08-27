@inline function _nu_ϵ⁺(f::Float64, g::Float64, ξ::Float64, lo, hi)::Float64
    if f < g
        return min(g - f, ξ*(f - lo))
    else
        return max(g - f, -ξ*(hi - f))
    end
end

@inline function _nu_ϵ⁻(f::Float64, g::Float64, ξ::Float64, lo, hi)::Float64
    if f < g
        return max(f - g, -ξ*(f - lo))
    else
        return min(f - g, ξ*(hi - f))
    end
end

_nuΦ⁺(α, f₋, f₀, f₊, d₋, d₀, d₊, ξ, lo, hi) = α*(f₀ + (d₀ - α)/(d₊ + d₀ + d₋)*(
                _nu_ϵ⁺(f₀, f₊, ξ, lo, hi)/(d₊ + d₀)*(d₀ + d₋ - α) +
                _nu_ϵ⁻(f₀, f₋, ξ, lo, hi)/(d₀ + d₋)*(d₊ + α)))

_nuΦ⁻(α, f₋, f₀, f₊, d₋, d₀, d₊, ξ, lo, hi) = α*(f₀ - (d₀ + α)/(d₊ + d₀ + d₋)*(
                _nu_ϵ⁺(f₀, f₊, ξ, lo, hi)/(d₊ + d₀)*(d₋ - α) +
                _nu_ϵ⁻(f₀, f₋, ξ, lo, hi)/(d₀ + d₋)*(d₀ + d₊ + α)))

"""
For `PFCNonUniform` the fourth argument is the displacement `vΔt`, a length,
not a Courant number: a non-uniform grid has no single Courant number to quote.
"""
function advect!(dest, src, p::PFCNonUniform, α, ws::PFCWorkspace)
    Δx, ξ, lo, hi = p.Δx, p.ξ, p.fmin, p.fmax
    n = length(Δx)
    acc = ws.accumulator
    copyto!(acc, src)

    if α > 0
        Φ = _nuΦ⁺(α, src[end], src[1], src[2], Δx[end], Δx[1], Δx[2], ξ[1], lo, hi)
        acc[1] = acc[1] - Φ/Δx[1]
        acc[2] = acc[2] + Φ/Δx[2]
        for i in 2:n-1
            Φ = _nuΦ⁺(α, src[i-1], src[i], src[i+1], Δx[i-1], Δx[i], Δx[i+1], ξ[i], lo, hi)
            acc[i] = acc[i] - Φ/Δx[i]
            acc[i+1] = acc[i+1] + Φ/Δx[i+1]
        end
        Φ = _nuΦ⁺(α, src[end-1], src[end], src[1], Δx[end-1], Δx[end], Δx[1], ξ[end], lo, hi)
        acc[end] = acc[end] - Φ/Δx[end]
        acc[1] = acc[1] + Φ/Δx[1]
    else
        Φ = _nuΦ⁻(α, src[end], src[1], src[2], Δx[end], Δx[1], Δx[2], ξ[1], lo, hi)
        acc[1] = acc[1] + Φ/Δx[1]
        acc[end] = acc[end] - Φ/Δx[end]
        for i in 2:n-1
            Φ = _nuΦ⁻(α, src[i-1], src[i], src[i+1], Δx[i-1], Δx[i], Δx[i+1], ξ[i], lo, hi)
            acc[i] = acc[i] + Φ/Δx[i]
            acc[i-1] = acc[i-1] - Φ/Δx[i-1]
        end
        Φ = _nuΦ⁻(α, src[end-1], src[end], src[1], Δx[end-1], Δx[end], Δx[1], ξ[end], lo, hi)
        acc[end] = acc[end] + Φ/Δx[end]
        acc[end-1] = acc[end-1] - Φ/Δx[end-1]
    end
    return copyto!(dest, acc)
end
