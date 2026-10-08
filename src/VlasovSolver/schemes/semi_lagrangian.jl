# The spline, its prefilter and its evaluation follow Interpolations.jl 0.16,
# `BSpline(Linear())` and `BSpline(Quadratic|Cubic(Periodic(OnCell())))` with
# `extrapolate(itp, Periodic(OnCell()))`, which this scheme called until 0.2:
# the same knots, the same weights in the same order of operations, and the same
# periodic image of a point. The prefilter solves the same system by a
# different factorisation, so the quadratic and cubic results moved in their
# last bits and the linear one not at all. `test/VlasovSolver/test_spline.jl`
# holds the two to each other.

function advect!(dest, src, s::SemiLagrangian, c, ws::SplineWorkspace)
    _validate(dest, src, s, c, ws)
    n = length(src)
    spline = s.spline
    coefficients = _prefilter!(ws, spline, src)
    c = float(c)
    # The cells whose point and stencil lie inside the knots go through a loop
    # with neither the periodic map nor a wrapped index, and the few at either
    # end through the general one. The arithmetic is the same in both.
    i₀, i₁ = _interior(c, n)
    _sample!(dest, coefficients, spline, c, n, 1:i₀-1)
    @inbounds for i = i₀:i₁
        dest[i] = _evaluate(coefficients, n, spline, i - c, Val(false))
    end
    _sample!(dest, coefficients, spline, c, n, i₁+1:n)
    return dest
end

function _sample!(dest, coefficients, spline, c, n, cells)
    u = _upper(spline, n)
    @inbounds for i in cells
        x = i - c
        # The old split between the interpolant and its periodic extrapolation:
        # only a point outside the knots is mapped.
        if !(1 ≤ x ≤ u)
            x = _periodic(x, _lower(spline, x), n)
        end
        dest[i] = _evaluate(coefficients, n, spline, x, Val(true))
    end
    return dest
end

# The cells `i₀:i₁` whose `x = i - c` surely lies in `[3, n - 2]`, where every
# spline's stencil is inside `1:n` and no point needs mapping: `i - c` is then
# above 3 and at most `n - 2` before rounding, which cannot cross either. Empty,
# as `(n + 1, n)`, when `c` takes every point outside.
function _interior(c, n)
    abs(c) < n || return (n + 1, n)
    fc = unsafe_trunc(Int, floor(c))
    i₀, i₁ = max(1, fc + 4), min(n, fc + n - 2)
    return i₀ ≤ i₁ ? (i₀, i₁) : (n + 1, n)
end

"""
    _prefilter!(ws, spline, src)

The B-spline coefficients of `src`, written into `ws.coefficients`.

The linear spline interpolates its data as they are, and gets the first point
again at `n + 1`, so that its last cell needs no wrapped index.

A quadratic or cubic spline interpolates the coefficients `c` that solve the
periodic system `M c = src`, `M` tridiagonal with `a` on the diagonal, `b` off
it and `b` in its corners `M[1,n]` and `M[n,1]` ([`_diagonals`](@ref)). With
`u = e₁ + eₙ`, `M = A + b·u·uᵀ`, where `A` is tridiagonal with `a − b` in its
two corners; `A` is strictly diagonally dominant, so needs no pivoting, and by
Sherman–Morrison

    c = y − b(y₁ + yₙ)·z,    y = A⁻¹src,    z = A⁻¹u/(1 + b·uᵀA⁻¹u).

The workspace holds `A`'s factorisation and `z`, so a step is a forward and a
backward sweep and one axpy, O(n), allocating nothing. Interpolations solves the
same system by a Woodbury update of an LU factorisation, rebuilt every step.
"""
function _prefilter!(ws::SplineWorkspace, ::LinearSpline, src)
    buf = ws.coefficients
    n = length(src)
    # Five-argument copyto! rather than a view, which on Julia 1.10 could be
    # heap-allocated depending on the host.
    copyto!(buf, 1, src, 1, n)
    buf[n+1] = src[1]
    return buf
end

function _prefilter!(ws::SplineWorkspace{T}, spline::Union{QuadraticSpline, CubicSpline},
                     src) where {T}
    c = ws.coefficients
    z = ws.z
    _, b = _diagonals(T, spline)
    copyto!(c, src)
    _thomas!(c, ws.rdiag, b)
    s = b*(c[1] + c[end])
    @inbounds @simd for i in eachindex(c, z)
        c[i] -= s*z[i]
    end
    return c
end

# The interval the knots span, `[_lower, _upper]`, which a point outside is
# mapped into. The linear spline's knots are `1, …, n + 1`, its last point a
# copy of the first; the others' are `1, …, n` and their cells reach half a
# cell beyond either end. Either way the period is `n`.
_upper(::LinearSpline, n) = n + 1
_upper(::AbstractSpline, n) = n
_lower(::LinearSpline, x) = one(x)
_lower(::AbstractSpline, x) = one(x)/2

# Interpolations' `periodic(x, l, u)`.
_periodic(x, l, n) = mod(x - l, oftype(x, n)) + l

# A knot index `1 - n ≤ j ≤ 2n`, taken back into `1:n`: Interpolations'
# `modrange`, without the division, for the range `_evaluate` asks for. Inside
# the interior, `_interior`, no index needs it.
_knot(j, n, ::Val{true}) = ifelse(j < 1, j + n, ifelse(j > n, j - n, j))
_knot(j, n, ::Val{false}) = j

"""
    _evaluate(coefficients, n, spline, x, wrap)

The spline with the given coefficients at `x`, which lies in
`[_lower, _upper]`, by Interpolations' `value_weights` and summed as its
`interp_getindex` sums, the last weight first. `wrap` is `Val(false)` only
where the stencil is known to lie inside `1:n`.
"""
@inline function _evaluate(buf, n, ::LinearSpline, x, wrap)
    f = floor(x)
    f = ifelse(x == n + 1, f - 1, f)    # x = n + 1 is in the last cell
    k = unsafe_trunc(Int, f)
    δ = x - f
    @inbounds return δ*buf[k+1] + (1 - δ)*buf[k]
end

@inline function _evaluate(c, n, ::QuadraticSpline, x, wrap)
    h = one(x)/2
    xh = x + h
    # Interpolations' `roundbounds`: a point on the upper edge goes to the cell
    # below it
    xm = ifelse(x < n + h, floor(xh), ceil(xh) - 1)
    k = unsafe_trunc(Int, xm)
    δ = x - xm
    w₁ = ((δ - h)*(δ - h))/2
    w₂ = (one(x)*3)/4 - δ*δ
    w₃ = ((δ + h)*(δ + h))/2
    @inbounds return w₃*c[_knot(k + 1, n, wrap)] + (w₂*c[_knot(k, n, wrap)] + w₁*c[_knot(k - 1, n, wrap)])
end

@inline function _evaluate(c, n, ::CubicSpline, x, wrap)
    f = floor(x)
    k = unsafe_trunc(Int, f)
    δ = x - f
    δ³ = δ*δ*δ
    v = 1 - δ
    v³ = v*v*v
    sixth = one(x)/6
    w₁ = sixth*v³
    w₂ = ((one(x)*2)/3 - δ*δ) + (one(x)/2)*δ³
    w₃ = ((one(x)*2)/3 - v*v) + (one(x)/2)*v³
    w₄ = sixth*δ³
    @inbounds return w₄*c[_knot(k + 2, n, wrap)] + (w₃*c[_knot(k + 1, n, wrap)] +
                     (w₂*c[_knot(k, n, wrap)] + w₁*c[_knot(k - 1, n, wrap)]))
end
