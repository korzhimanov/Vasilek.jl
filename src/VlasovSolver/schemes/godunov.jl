# The limiter's argument, the upwind difference over the downwind one: 0 where
# the upwind side is flat, 2 where only the downwind side is.
#
# `iszero` rather than `≈ 0.0`, which with its default tolerances is the same
# test for every Float64, NaN and ±0 included, but gets there through
# `isfinite` checks and a NaN-aware `max`: in the vectorized loop, twice what
# the rest of the step costs. And `ifelse` rather than `if`, so that the loop
# vectorizing does not rest on LLVM turning the branches into selects. The
# division is evaluated whether or not it is used.
@inline function _ratio(f, i, i⁻, i⁺)
    Δ⁻ = f[i] - f[i⁻]
    Δ⁺ = f[i⁺] - f[i]
    return ifelse(iszero(Δ⁻), 0.0, ifelse(iszero(Δ⁺), 2.0, Δ⁻/Δ⁺))
end

# The flux through the face between cell i and its upwind neighbour i⁻, which
# leaves i⁻ and enters i: the linear reconstruction in i⁻, averaged over the
# strip that crosses the face in one step, whose midpoint is |c|Δx/2 upwind of
# the face, (1 - |c|)Δx/2 downwind of the cell's centre. `@inline`, like
# `_ratio`: left to itself, the inliner declined `_Φ` with a limiter on Julia
# 1.13, and `_ratio` with or without one on 1.10, where the call stayed even when
# `NoLimiter` discards its result, because its bounds checks can throw. Either
# way the loop made two calls per cell, vectorized nothing, and bounds-checked
# every read inside them.
@inline function _Φ(f, i, i⁻, i⁻², c, ::PiecewiseLinear, limiter)
    a = abs(c)
    return a*(f[i⁻] + limiter(_ratio(f, i⁻, i⁻², i))*0.5*(1 - a)*(f[i] - f[i⁻]))
end

# One `advect!` for every reconstruction, so the step is validated in one place;
# the sweep dispatches on the reconstruction below it. A method per
# reconstruction would be a `_validate` call per reconstruction, each to be kept
# in step with that function's signature -- and a signature change made
# elsewhere merges without a conflict and leaves a call it missed failing at run
# time.
function advect!(dest, src, g::Godunov, c, ws)
    _validate(dest, src, g, c, ws)
    return _godunov!(dest, src, c, g.reconstruction, g.limiter)
end

# Each face's flux is evaluated once: the loop computes the flux through face
# i+½, uses it in cell i, and carries it to cell i+1, where it is the flux
# through face i-½. Evaluating both faces of every cell, as this used to, did
# the limiter's two divisions twice. `Φ₋` and `Φ₊` are the fluxes through faces
# i-½ and i+½, and `Φ½` the one through the periodic face between cells n and 1.
# A cell's update is (value + inflow) - outflow, the expression the per-cell
# form evaluated, so the result is the same to the bit.
#
# No `@simd`: LLVM vectorizes the carried flux as a first-order recurrence
# without it, and the kernel has no reduction for it to reassociate.
function _godunov!(dest, src, c, r, l)
    n = length(dest)
    if c > 0
        # face i+½'s upwind cell is i
        Φ½ = _Φ(src, 1, n, n-1, c, r, l)
        Φ₋ = _Φ(src, 2, 1, n, c, r, l)
        dest[1] = src[1] + Φ½ - Φ₋
        @inbounds for i = 2:n-1
            Φ₊ = _Φ(src, i+1, i, i-1, c, r, l)
            dest[i] = src[i] + Φ₋ - Φ₊
            Φ₋ = Φ₊
        end
        dest[n] = src[n] + Φ₋ - Φ½
    else
        # the mirror image: face i+½'s upwind cell is i+1
        Φ½ = _Φ(src, n, 1, 2, c, r, l)
        Φ₋ = Φ½
        @inbounds for i = 1:n-2
            Φ₊ = _Φ(src, i, i+1, i+2, c, r, l)
            dest[i] = src[i] + Φ₊ - Φ₋
            Φ₋ = Φ₊
        end
        Φ₊ = _Φ(src, n-1, n, 1, c, r, l)
        dest[n-1] = src[n-1] + Φ₊ - Φ₋
        dest[n] = src[n] + Φ½ - Φ₊
    end
    return dest
end

# The constant reconstruction's flux, |c| times the upwind cell, is one
# multiplication. Evaluating it for both faces of every cell costs less than
# carrying it, which takes two vector shuffles per four cells: the loop above
# measured 8% to 27% slower on this reconstruction, on Julia 1.10 and 1.13.
function _godunov!(dest, src, c, ::PiecewiseConstant, _)
    n = length(dest)
    a = abs(c)
    if c > 0
        dest[1] = src[1] + a*src[n] - a*src[1]
        @inbounds @simd for i = 2:n
            dest[i] = src[i] + a*src[i-1] - a*src[i]
        end
    else
        @inbounds @simd for i = 1:n-1
            dest[i] = src[i] + a*src[i+1] - a*src[i]
        end
        dest[n] = src[n] + a*src[1] - a*src[n]
    end
    return dest
end
