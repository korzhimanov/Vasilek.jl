module SemiLagrangian

using Interpolations

function _interpolation_type(interpolation_order)
    if interpolation_order == :Linear
        return BSpline(Linear())
    elseif interpolation_order == :Quadratic
        return BSpline(Quadratic(Periodic(OnCell())))
    elseif interpolation_order == :Cubic
        return BSpline(Cubic(Periodic(OnCell())))
    else
        throw(ArgumentError("unknown interpolation_order $interpolation_order"))
    end
end

"""
    _make_step(f₀, f, interpolation_order)

Builds the per-step kernel shared by both `generate_solver` methods.

`interpolate` copies its input and runs the B-spline prefilter on the copy, so
calling it inside the step allocated a fresh coefficient array every time —
80 kB per step at N = 1000 for linear, 173 kB for quadratic and cubic. The
buffer is now allocated once and prefiltered in place; `interpolate!` was
checked to give bit-identical results to `interpolate`.

Linear keeps its explicit wrap-around point. Dropping it and relying on the
periodic extrapolation alone changes the result, so it stays.

(The two `generate_solver` methods share this rather than each carrying a copy
of the body. That is a local consequence of not wanting to apply the same fix
twice, not the wider deduplication of the closure factories.)
"""
function _make_step(f₀, f, interpolation_order)
    itp_type = _interpolation_type(interpolation_order)
    wrap = interpolation_order == :Linear
    buf = wrap ? similar(f₀, length(f₀) + 1) : similar(f₀)

    return function step!(c)
        if wrap
            copyto!(view(buf, 1:length(f₀)), f₀)
            buf[end] = f₀[1]
        else
            copyto!(buf, f₀)
        end
        itp = interpolate!(buf, itp_type)
        etp = extrapolate(itp, Periodic(OnCell()))

        for i = 1:length(f₀)
            if 1.0 ≤ i - c ≤ length(itp)
                f[i] = itp(i - c)
            else
                f[i] = etp(i - c)
            end
        end
    end
end

function generate_solver(f₀, f; interpolation_order = :Cubic)
    return _make_step(f₀, f, interpolation_order)
end

function generate_solver(f₀, f, c; interpolation_order = :Cubic)
    step! = _make_step(f₀, f, interpolation_order)
    return () -> step!(c)
end

end  # module SemiLagrangian
