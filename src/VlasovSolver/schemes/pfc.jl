@inline function _ϵ⁺(f, g, fmin, fmax)
    if f < g
        return min(g - f, 2*(f - fmin))
    else
        return max(g - f, -2*(fmax - f))
    end
end

@inline function _ϵ⁻(f, g, fmin, fmax)
    if f < g
        return max(f - g, -2*(f - fmin))
    else
        return min(f - g, 2*(fmax - f))
    end
end

_Φ⁺(f, i, i⁻, i⁺, c, lo, hi) = c*(f[i] + (1 - c)/3*(
                _ϵ⁺(f[i], f[i⁺], lo, hi)/2*(2 - c) +
                _ϵ⁻(f[i], f[i⁻], lo, hi)/2*(1 + c)))

_Φ⁻(f, i, i⁻, i⁺, c, lo, hi) = c*(f[i] - (1 + c)/3*(
                _ϵ⁺(f[i], f[i⁺], lo, hi)/2*(1 - c) +
                _ϵ⁻(f[i], f[i⁻], lo, hi)/2*(2 + c)))

function advect!(dest, src, p::PFC{T,Checked}, c, ws) where {T,Checked}
    _validate(dest, src, p, c, ws)
    n = length(dest)
    lo, hi = p.fmin, p.fmax
    if Checked
        # A minimum/maximum pass per call, about 7% of the step at N = 10000 on
        # Julia 1.13 and 17% on 1.10. Compiled away entirely when the scheme is
        # built with checked = false.
        _check_bounds(src, lo, hi)
    end

    # Each cell takes the difference of the fluxes through its two faces,
    # `src[i] + (Φleft - Φright)` (signed rightwards), as `PFCNonUniform` does,
    # so that equal fluxes cancel before they reach `f`. On a constant line
    # every face carries the same `c*f` to the bit, and the line stays constant
    # to the bit. The update was `(src[i] + Φin) - Φout`, which rounds twice:
    # from a constant line at `fmax = 0.2631…` it left a cell one ulp above
    # `fmax` at 320 of 2001 Courant numbers across [-1, 1], and the next
    # checked call refused it.
    if c > 0
        # `_Φ⁺` at cell i is the flux through its right face.
        Φ₁ = _Φ⁺(src, 1, n, 2, c, lo, hi)
        Φ = Φ₁
        for i in 2:n-1
            Φright = _Φ⁺(src, i, i-1, i+1, c, lo, hi)
            dest[i] = src[i] + (Φ - Φright)
            Φ = Φright
        end
        Φₙ = _Φ⁺(src, n, n-1, 1, c, lo, hi)
        dest[n] = src[n] + (Φ - Φₙ)
        dest[1] = src[1] + (Φₙ - Φ₁)
    else
        # `_Φ⁻` at cell i is the flux through its left face.
        Φ₁ = _Φ⁻(src, 1, n, 2, c, lo, hi)
        Φ = Φ₁
        for i in 2:n-1
            Φright = _Φ⁻(src, i, i-1, i+1, c, lo, hi)
            dest[i-1] = src[i-1] + (Φ - Φright)
            Φ = Φright
        end
        Φₙ = _Φ⁻(src, n, n-1, 1, c, lo, hi)
        dest[n-1] = src[n-1] + (Φ - Φₙ)
        dest[n] = src[n] + (Φₙ - Φ₁)
    end
    return dest
end
