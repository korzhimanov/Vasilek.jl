function _ratio(f, i, i⁻, i⁺)
    if f[i] - f[i⁻] ≈ 0.0
        return 0.0
    elseif f[i⁺] - f[i] ≈ 0.0
        return 2.0
    else
        return (f[i] - f[i⁻])/(f[i⁺] - f[i])
    end
end

_Φ(f, i, i⁻, i⁻², c, ::PiecewiseConstant, limiter) = abs(c)*f[i⁻]

function _Φ(f, i, i⁻, i⁻², c, ::PiecewiseLinear, limiter)
    return abs(c)*(f[i⁻] + limiter(_ratio(f, i⁻, i⁻², i))*0.5*(f[i] - f[i⁻]))
end

@inline function _godunov(f, i, i⁻, i⁺, i⁻², c, r, l)
    return f[i] + _Φ(f, i, i⁻, i⁻², c, r, l) - _Φ(f, i⁺, i, i⁻, c, r, l)
end

function advect!(dest, src, g::Godunov, c, ws)
    _validate(dest, src, g, ws)
    n = length(dest)
    r, l = g.reconstruction, g.limiter
    if c > 0
        dest[1] = _godunov(src, 1, n, 2, n-1, c, r, l)
        dest[2] = _godunov(src, 2, 1, 3, n, c, r, l)
        @inbounds @simd for i = 3:n-1
            dest[i] = _godunov(src, i, i-1, i+1, i-2, c, r, l)
        end
        dest[n] = _godunov(src, n, n-1, 1, n-2, c, r, l)
    else
        # the negative direction is the positive one with the neighbours swapped
        dest[1] = _godunov(src, 1, 2, n, 3, c, r, l)
        @inbounds @simd for i = 2:n-2
            dest[i] = _godunov(src, i, i+1, i-1, i+2, c, r, l)
        end
        dest[n-1] = _godunov(src, n-1, n, n-2, 1, c, r, l)
        dest[n] = _godunov(src, n, 1, n-1, 2, c, r, l)
    end
    return dest
end
