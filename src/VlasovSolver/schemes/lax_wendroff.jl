@inline function _lax_wendroff(f, i, i⁻, i⁺, c)
    return f[i] - 0.5*c*(f[i⁺] - f[i⁻]) + 0.5*c^2*(f[i⁺] - 2*f[i] + f[i⁻])
end

function advect!(dest, src, scheme::LaxWendroff, c, ws)
    _validate(dest, src, scheme, ws)
    n = length(dest)
    dest[1] = _lax_wendroff(src, 1, n, 2, c)
    @inbounds @simd for i = 2:n-1
        dest[i] = _lax_wendroff(src, i, i-1, i+1, c)
    end
    dest[n] = _lax_wendroff(src, n, n-1, 1, c)
    return dest
end
