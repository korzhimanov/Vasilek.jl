function advect!(dest, src, s::SemiLagrangian{LinearSpline}, c, ws::SplineWorkspace)
    n = length(src)
    buf = ws.buffer
    # explicit wrap-around point: relying on the periodic extrapolation alone
    # changes the result. Five-argument copyto! rather than a view, which on
    # Julia 1.10 could be heap-allocated depending on the host.
    copyto!(buf, 1, src, 1, n)
    buf[end] = src[1]
    return _sample!(dest, n, buf, _bspline(s.spline), c)
end

function advect!(dest, src, s::SemiLagrangian, c, ws::SplineWorkspace)
    n = length(src)
    buf = ws.buffer
    copyto!(buf, src)
    return _sample!(dest, n, buf, _bspline(s.spline), c)
end

function _sample!(dest, n, buf, itp_type, c)
    itp = interpolate!(buf, itp_type)
    etp = extrapolate(itp, Periodic(OnCell()))
    for i = 1:n
        if 1.0 ≤ i - c ≤ length(itp)
            dest[i] = itp(i - c)
        else
            dest[i] = etp(i - c)
        end
    end
    return dest
end
