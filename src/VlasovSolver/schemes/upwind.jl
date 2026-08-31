@inline _upwind⁺(f, i, i⁻, c) = f[i] - c*(f[i] - f[i⁻])
@inline _upwind⁻(f, i, i⁺, c) = f[i] - c*(f[i⁺] - f[i])

function advect!(dest, src, scheme::Upwind, c, ws)
    _validate(dest, src, scheme, ws)
    n = length(dest)
    if c > 0
        dest[1] = _upwind⁺(src, 1, n, c)
        @inbounds @simd for i = 2:n
            dest[i] = _upwind⁺(src, i, i-1, c)
        end
    else
        @inbounds @simd for i = 1:n-1
            dest[i] = _upwind⁻(src, i, i+1, c)
        end
        dest[n] = _upwind⁻(src, n, 1, c)
    end
    return dest
end
