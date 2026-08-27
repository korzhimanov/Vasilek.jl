module PFCNonUniform
export make_advect_1D!

"""
    slope_limit(r)

PFC slope-limiter coefficient for a cell triple whose smallest-to-largest
spacing ratio is `r ∈ (0, 1]`. Equals 2 on a locally uniform grid.
"""
slope_limit(r) = (1.0 + r)*(1.0 + 2r)/(3.0 + (r - 1.0/r)^2)

"""
    make_advect_1D!(Δx; fₘᵢₙ, fₘₐₓ)

Positive Flux Conservative advection on a static non-uniform 1D grid with
periodic boundaries. Returns `advect_1D!(f, α)`.

`fₘᵢₙ` and `fₘₐₓ` are required, as in [`PFC`](@ref): the limiter clamps
against them, and by Liouville's theorem the caller must supply the global
bounds of the *initial* condition. They used to be hardcoded to 0 and 1 --
`ϵ⁺` read `ξ*(1 - f)` -- which is silently wrong for any `f` exceeding 1.

The limiter coefficient is computed **per cell triple** rather than once from
the global `minimum(Δx)/maximum(Δx)`. With a global value, a single refined
region tightens the limiter across the whole domain: on a grid whose spacing
alternates between 0.1 and 0.2 the global ratio gives ξ ≈ 0.571 everywhere,
against ξ = 2 on a uniform grid — a limiter three and a half times more
restrictive even in the regions that are locally uniform. On a uniform grid
the two agree exactly, so this changes nothing there.
"""
function make_advect_1D!(Δx_::Vector{Float64}; fₘᵢₙ, fₘₐₓ)
    Δx = copy(Δx_)
    n = length(Δx)

    # ξ depends only on the grid, so it is computed once, per cell triple,
    # with periodic wrap.
    ξ = similar(Δx)
    for i in eachindex(Δx)
        d₋ = Δx[i == 1 ? n : i-1]
        d₊ = Δx[i == n ? 1 : i+1]
        ξ[i] = slope_limit(min(d₋, Δx[i], d₊)/max(d₋, Δx[i], d₊))
    end

    f_tmp = similar(Δx)

    function ϵ⁺(f::Float64, g::Float64, ξ::Float64)::Float64
        if f < g
            return min(g - f, ξ*(f - fₘᵢₙ))
        else
            return max(g - f, -ξ*(fₘₐₓ - f))
        end
    end

    function ϵ⁻(f::Float64, g::Float64, ξ::Float64)::Float64
        if f < g
            return max(f - g, -ξ*(f - fₘᵢₙ))
        else
            return min(f - g, ξ*(fₘₐₓ - f))
        end
    end

    # scalar arguments: the previous signature took three-element slices, so
    # every cell of every step allocated two temporary vectors
    Φ⁺(α, f₋, f₀, f₊, d₋, d₀, d₊, ξ) = α*(f₀ + (d₀ - α)/(d₊ + d₀ + d₋)*(
                    ϵ⁺(f₀, f₊, ξ)/(d₊ + d₀)*(d₀ + d₋ - α) +
                    ϵ⁻(f₀, f₋, ξ)/(d₀ + d₋)*(d₊ + α)))

    Φ⁻(α, f₋, f₀, f₊, d₋, d₀, d₊, ξ) = α*(f₀ - (d₀ + α)/(d₊ + d₀ + d₋)*(
                    ϵ⁺(f₀, f₊, ξ)/(d₊ + d₀)*(d₋ - α) +
                    ϵ⁻(f₀, f₋, ξ)/(d₀ + d₋)*(d₀ + d₊ + α)))

    function advect_1D!(f, α::Float64)
        copyto!(f_tmp, f)
        if α > 0
            Φ = Φ⁺(α, f[end], f[1], f[2], Δx[end], Δx[1], Δx[2], ξ[1])
            f_tmp[1] = f_tmp[1] - Φ/Δx[1]
            f_tmp[2] = f_tmp[2] + Φ/Δx[2]
            for i in 2:n-1
                Φ = Φ⁺(α, f[i-1], f[i], f[i+1], Δx[i-1], Δx[i], Δx[i+1], ξ[i])
                f_tmp[i] = f_tmp[i] - Φ/Δx[i]
                f_tmp[i+1] = f_tmp[i+1] + Φ/Δx[i+1]
            end
            Φ = Φ⁺(α, f[end-1], f[end], f[1], Δx[end-1], Δx[end], Δx[1], ξ[end])
            f_tmp[end] = f_tmp[end] - Φ/Δx[end]
            f_tmp[1] = f_tmp[1] + Φ/Δx[1]
        else
            Φ = Φ⁻(α, f[end], f[1], f[2], Δx[end], Δx[1], Δx[2], ξ[1])
            f_tmp[1] = f_tmp[1] + Φ/Δx[1]
            f_tmp[end] = f_tmp[end] - Φ/Δx[end]
            for i in 2:n-1
                Φ = Φ⁻(α, f[i-1], f[i], f[i+1], Δx[i-1], Δx[i], Δx[i+1], ξ[i])
                f_tmp[i] = f_tmp[i] + Φ/Δx[i]
                f_tmp[i-1] = f_tmp[i-1] - Φ/Δx[i-1]
            end
            Φ = Φ⁻(α, f[end-1], f[end], f[1], Δx[end-1], Δx[end], Δx[1], ξ[end])
            f_tmp[end] = f_tmp[end] + Φ/Δx[end]
            f_tmp[end-1] = f_tmp[end-1] - Φ/Δx[end-1]
        end
        copyto!(f, f_tmp)
    end

    return advect_1D!
end

end
