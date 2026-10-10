# Factored out deliberately. This used to be written twice inside the solver,
# and the second copy branched on the outer loop index while claiming to
# differentiate at the inner one, so the two "different" derivatives were always
# equal and the collision operator computed something else entirely.
"""
    ∂f∂v(f, v, k)

Derivative of `f` with respect to `v` at index `k`: centred in the interior,
one-sided at either end.
"""
function ∂f∂v(f, v, k)
    if k == 1
        return (f[2] - f[1])/(v[2] - v[1])
    elseif k == length(f)
        return (f[end] - f[end-1])/(v[end] - v[end-1])
    else
        return (f[k+1] - f[k-1])/(v[k+1] - v[k-1])
    end
end

# Where the model comes from. The bracket is the `x` row of the Landau tensor
# `(|u|²I − uu)/|u|³` with `∂f/∂v⊥ = −(v⊥/Tₜ) f`: its `xx` component gives the
# first two terms, its `xy` and `xz` components the drag `−f f′ (v − v′)/Tₜ`.
# Both carry `u⊥²/|u|³`, which `Φ` closes with `u⊥² ≈ 2Tₜ`, the operator's
# long-standing estimate; for `|u| ≫ √(2Tₜ)` it gives back the former kernel
# `2Tₜ/|u|³`. Averaging over both particles' transverse Maxwellians would give
# `4Tₜ` instead, but the choice sets only the kernel's width: for a Maxwellian
# at `T` the bracket is `f f′ (v − v′)(1/T − 1/Tₜ)` under any closure, so the
# equilibrium is the Maxwellian at `Tₜ` regardless.
"""
    Landau1P(A; L = 20.0, Tₜ = 1.0)

The Landau collision operator in one velocity dimension, the two transverse
ones a Maxwellian bath at temperature `Tₜ`:

    ∂f/∂t = −∂F/∂v,  F(v) = L·A ∫ Φ(v − v′) [f ∂f′ − f′ ∂f − f f′ (v − v′)/Tₜ] dv′,
    Φ(u) = 2Tₜ/(u² + 2Tₜ)^(3/2),

where `A = 4πe⁴N₀/(m²v₀³ω)` and `L` is the Coulomb logarithm. `Φ` closes the
transverse `u⊥²` with `2Tₜ`; both particles' Maxwellians would give `4Tₜ`.
That sets the kernel's width `√(2Tₜ)`, which the grid must resolve, but not
the equilibrium: the Maxwellian at `Tₜ`, by default the thermal one, is
stationary, any other temperature
relaxes to it, density and momentum are conserved, and the entropy relative to
that Maxwellian does not increase.

On the grid `F` is taken at the half-points from `f½ = (fᵢ + fᵢ₊₁)/2`,
`f′½ = (fᵢ₊₁ − fᵢ)/(vᵢ₊₁ − vᵢ)` and the nodal `fⱼ`, `f′ⱼ` weighted by the cell
widths `wⱼ`. It is zero at the window's ends, and
`destᵢ = srcᵢ − Δt (Fᵢ₊½ − Fᵢ₋½)/wᵢ`. So `Σ w f` is conserved to round-off;
momentum to second order in the spacing, and exactly on a uniform grid while
`f` vanishes at the ends; the Maxwellian at `Tₜ` is stationary to second
order. Forward Euler is stable for `Δt ≲ Δv²/(4·L·A·max f)`. A step is O(n²).

Not exported: reach it as `Vasilek.Collisions.Landau1P`.
"""
struct Landau1P{T<:AbstractFloat} <: AbstractCollisionOperator
    A::T
    L::T
    Tₜ::T
end
Landau1P(A; L = 20.0, Tₜ = 1.0) = Landau1P(promote(float(A), float(L), float(Tₜ))...)

# `F[k]` is the flux through the half-point between nodes k − 1 and k, so `F[1]`
# and `F[n+1]` are the walls at the window's ends; `df` is the nodal derivative,
# taken once per call.
struct Landau1PWorkspace{T}
    F::Vector{T}
    df::Vector{T}
end
workspace(::Landau1P{T}, n::Integer, ::Type{S} = T) where {T, S} =
    Landau1PWorkspace(Vector{S}(undef, n + 1), Vector{S}(undef, n))

function collide!(dest, src, op::Landau1P, v, Δt, ws::Landau1PWorkspace)
    n = length(src)
    length(dest) == n == length(v) || throw(DimensionMismatch(
        "collide! needs dest, src and v of one length, got $(length(dest)), " *
        "$n and $(length(v))"))
    F, df = ws.F, ws.df
    # Both are written while `src` is still being read.
    (Base.mightalias(src, F) || Base.mightalias(src, df)) && _err_landau_alias()
    n < 2 && return copyto!(dest, src)       # no half-point, nothing can flow
    # One type for the sum, the widest of the scratch's, the line's and the
    # grid's: an accumulator that widened on its first `+=` would be a Union.
    R = promote_type(eltype(F), float(eltype(src)), float(eltype(v)))
    Tₜ = R(op.Tₜ)
    σ² = 2Tₜ
    κ = inv(Tₜ)
    LAσ² = R(op.L*op.A)*σ²                   # Φ's numerator, out of the sum
    for j in 1:n
        df[j] = ∂f∂v(src, v, j)
    end
    F[1] = F[n+1] = zero(eltype(F))
    for i in 1:n-1
        v½ = (v[i] + v[i+1])/2
        f½ = (src[i] + src[i+1])/2
        f′½ = (src[i+1] - src[i])/(v[i+1] - v[i])
        s = zero(R)
        for j in 1:n
            u = v½ - v[j]
            q = u*u + σ²
            # wⱼ Φ(u)/σ² times the bracket, whose drag f½ fⱼ u/Tₜ rides with
            # fⱼ f′½; (u² + σ²)^(3/2) as q·√q, a square root where `^1.5`
            # would be a call to `pow`
            s += _width(v, j)*(f½*df[j] - src[j]*(f′½ + f½*u*κ))/(q*sqrt(q))
        end
        F[i+1] = LAσ²*s
    end
    for i in 1:n
        dest[i] = src[i] - Δt*(F[i+1] - F[i])/_width(v, i)
    end
    return dest
end

@noinline _err_landau_alias() = throw(ArgumentError(
    "collide! requires src not to share memory with the workspace: Landau1P " *
    "writes the derivative and the fluxes there before it is done reading src"))
