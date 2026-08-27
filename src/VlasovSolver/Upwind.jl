module Upwind

@inline _flux⁺(f, i, i⁻, c) = f[i] - c*(f[i] - f[i⁻])
@inline _flux⁻(f, i, i⁺, c) = f[i] - c*(f[i⁺] - f[i])

"""
    _advect!(dest, src, c)

First-order upwind advection with periodic boundaries, one step.

The single numerical kernel. Both `generate_solver` methods are wrappers over
it: the only thing that ever differed between them was whether the Courant
number was captured at construction or passed at call time, and the kernel was
written out twice to express that.
"""
function _advect!(dest, src, c)
    n = length(dest)
    if c > 0
        dest[1] = _flux⁺(src, 1, n, c)
        @inbounds @simd for i = 2:n
            dest[i] = _flux⁺(src, i, i-1, c)
        end
    else
        @inbounds @simd for i = 1:n-1
            dest[i] = _flux⁻(src, i, i+1, c)
        end
        dest[n] = _flux⁻(src, n, 1, c)
    end
    return dest
end

generate_solver(f₀, f) = c -> _advect!(f, f₀, c)
generate_solver(f₀, f, c) = () -> _advect!(f, f₀, c)

end # module
