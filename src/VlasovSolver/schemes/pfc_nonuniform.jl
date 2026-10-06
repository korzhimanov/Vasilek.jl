@inline function _nu_ϵ⁺(f, g, ξ, lo, hi)
    if f < g
        return min(g - f, ξ*(f - lo))
    else
        return max(g - f, -ξ*(hi - f))
    end
end

@inline function _nu_ϵ⁻(f, g, ξ, lo, hi)
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
function advect!(dest, src, p::PFCNonUniform{T,Checked}, α,
                 ws::PFCWorkspace) where {T,Checked}
    _validate(dest, src, p, α, ws)
    Δx, ξ, lo, hi = p.Δx, p.ξ, p.fmin, p.fmax
    if Checked
        # PFC's minimum/maximum pass, about 5% of this step at N = 10000 on
        # Julia 1.13 and 10% on 1.10. Compiled away when checked = false.
        _check_bounds(src, lo, hi)
    end
    n = length(Δx)
    acc = ws.accumulator

    # Each cell takes the difference of the fluxes through its two faces,
    # `src[i] + (Φleft - Φright)/Δx[i]` (signed rightwards), so that equal
    # fluxes cancel before they reach `f`. On a constant line every face
    # carries the same `α*f` to the bit, and the line stays constant to the
    # bit. The update was `(src[i] + Φin/Δx[i]) - Φout/Δx[i]`, which rounds
    # twice: from a constant line at `fmax` it left a cell one ulp above
    # `fmax` in 691 of the 2406 cases `test_nonuniform_advection.jl` runs,
    # uniform grids among them, and the next checked call refused it.
    if α > 0
        # `_nuΦ⁺` at cell i is the flux through its right face.
        Φ₁ = _nuΦ⁺(α, src[end], src[1], src[2], Δx[end], Δx[1], Δx[2], ξ[1], lo, hi)
        Φ = Φ₁
        for i in 2:n-1
            Φright = _nuΦ⁺(α, src[i-1], src[i], src[i+1], Δx[i-1], Δx[i], Δx[i+1], ξ[i], lo, hi)
            acc[i] = src[i] + (Φ - Φright)/Δx[i]
            Φ = Φright
        end
        Φₙ = _nuΦ⁺(α, src[end-1], src[end], src[1], Δx[end-1], Δx[end], Δx[1], ξ[end], lo, hi)
        acc[end] = src[end] + (Φ - Φₙ)/Δx[end]
        acc[1] = src[1] + (Φₙ - Φ₁)/Δx[1]
    else
        # `_nuΦ⁻` at cell i is the flux through its left face.
        Φ₁ = _nuΦ⁻(α, src[end], src[1], src[2], Δx[end], Δx[1], Δx[2], ξ[1], lo, hi)
        Φ = Φ₁
        for i in 2:n-1
            Φright = _nuΦ⁻(α, src[i-1], src[i], src[i+1], Δx[i-1], Δx[i], Δx[i+1], ξ[i], lo, hi)
            acc[i-1] = src[i-1] + (Φ - Φright)/Δx[i-1]
            Φ = Φright
        end
        Φₙ = _nuΦ⁻(α, src[end-1], src[end], src[1], Δx[end-1], Δx[end], Δx[1], ξ[end], lo, hi)
        acc[end-1] = src[end-1] + (Φ - Φₙ)/Δx[end-1]
        acc[end] = src[end] + (Φₙ - Φ₁)/Δx[end]
    end
    return copyto!(dest, acc)
end
