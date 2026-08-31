@inline function _ϵ⁺(f, g, fmin, fmax)
    if f < g
        return min(g - f, 2.0*(f - fmin))
    else
        return max(g - f, -2.0*(fmax - f))
    end
end

@inline function _ϵ⁻(f, g, fmin, fmax)
    if f < g
        return max(f - g, -2.0*(f - fmin))
    else
        return min(f - g, 2.0*(fmax - f))
    end
end

_Φ⁺(f, i, i⁻, i⁺, c, lo, hi) = c*(f[i] + (1.0 - c)/3.0*(
                _ϵ⁺(f[i], f[i⁺], lo, hi)/2.0*(2.0 - c) +
                _ϵ⁻(f[i], f[i⁻], lo, hi)/2.0*(1.0 + c)))

_Φ⁻(f, i, i⁻, i⁺, c, lo, hi) = c*(f[i] - (1.0 + c)/3.0*(
                _ϵ⁺(f[i], f[i⁺], lo, hi)/2.0*(1.0 - c) +
                _ϵ⁻(f[i], f[i⁻], lo, hi)/2.0*(2.0 + c)))

function advect!(dest, src, p::PFC{T,Checked}, c, ws) where {T,Checked}
    _validate(dest, src, p, ws)
    n = length(dest)
    lo, hi = p.fmin, p.fmax
    if Checked
        # A minimum/maximum pass per call, about 13% of the step at N = 10000.
        # Compiled away entirely when the scheme is built with checked = false.
        @assert lo ≤ minimum(src) "fmin = $lo exceeds minimum(src) = $(minimum(src))"
        @assert maximum(src) ≤ hi "fmax = $hi is below maximum(src) = $(maximum(src))"
    end

    if c > 0
        Φ = _Φ⁺(src, 1, n, 2, c, lo, hi)
        dest[1] = src[1] - Φ
        dest[2] = src[2] + Φ
        for i in 2:n-1
            Φ = _Φ⁺(src, i, i-1, i+1, c, lo, hi)
            dest[i] = dest[i] - Φ
            dest[i+1] = src[i+1] + Φ
        end
        Φ = _Φ⁺(src, n, n-1, 1, c, lo, hi)
        dest[n] = dest[n] - Φ
        dest[1] = dest[1] + Φ
    else
        Φ = _Φ⁻(src, 1, n, 2, c, lo, hi)
        dest[1] = src[1] + Φ
        dest[n] = src[n] - Φ
        for i in 2:n-1
            Φ = _Φ⁻(src, i, i-1, i+1, c, lo, hi)
            dest[i] = src[i] + Φ
            dest[i-1] = dest[i-1] - Φ
        end
        Φ = _Φ⁻(src, n, n-1, 1, c, lo, hi)
        dest[n] = dest[n] + Φ
        dest[n-1] = dest[n-1] - Φ
    end
    return dest
end
