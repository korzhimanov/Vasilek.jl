module LaxWendroff

@inline function _flux(f, i, i⁻, i⁺, c)
    return f[i] - 0.5*c*(f[i⁺] - f[i⁻]) + 0.5*c^2*(f[i⁺] - 2*f[i] + f[i⁻])
end

"""
    _advect!(dest, src, c)

Lax–Wendroff advection with periodic boundaries, one step.

The single numerical kernel; both `generate_solver` methods wrap it.

Sizing from `length(dest)` also fixes a latent bug: the `(h, h₀)` methods used
to take their bounds from `length(f)`, the array captured at construction,
so passing a differently sized pair read and wrote out of range. The same
mistake was present in `PFC` and is fixed there too. No test exercised it.
"""
function _advect!(dest, src, c)
    n = length(dest)
    dest[1] = _flux(src, 1, n, 2, c)
    @inbounds @simd for i = 2:n-1
        dest[i] = _flux(src, i, i-1, i+1, c)
    end
    dest[n] = _flux(src, n, n-1, 1, c)
    return dest
end

function generate_solver(f₀, f)
    solve!(c) = _advect!(f, f₀, c)
    solve!(h, h₀, c) = _advect!(h, h₀, c)
    return solve!
end

function generate_solver(f₀, f, c)
    solve!() = _advect!(f, f₀, c)
    solve!(h, h₀) = _advect!(h, h₀, c)
    return solve!
end

end # module
