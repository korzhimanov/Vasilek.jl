module PoissonFourier1D

using FFTW, LinearAlgebra
import ..workspace

export PoissonFFT1D, solve!

"""
    PoissonFFT1D(n, Δx; derivative = :centered)

Spectral Poisson solve on a periodic grid of `n` points spaced `Δx`. `solve!`
returns the field `e` with `∂e/∂x = ρ`: in the units of
`docs/src/normalization.md`, pass `ρ = nₑ − nᵢ`. The zero mode of `ρ` is
discarded, so a net charge is ignored rather than producing a divergent field.

`derivative` chooses how `e` is taken from the potential:

  * `:centered` (the default, and what every quoted number in the repository
    was measured with): solve for `φ` spectrally, then a centred difference.
    A mode's field comes out scaled by `sin(kΔx)/(kΔx)`.
  * `:spectral`: `ê = −iρ̂/k` directly, exact for every resolved mode, one
    loop shorter. The Nyquist mode of an even `n` has no well-defined
    derivative and is set to zero, as the centred difference also gives.

The value holds no buffers and can be shared between tasks; the scratch
belongs to a [`workspace`](@ref), one per task.
"""
struct PoissonFFT1D{T<:AbstractFloat, C<:AbstractVector}
    n::Int
    Δx::T
    # multiplies ρ̂ mode by mode: real -1/k² for the centred form, which keeps
    # it bit-identical to the 0.1 solver, and complex -i/k for the spectral one
    coefficient::C
    spectral::Bool
end

function PoissonFFT1D(n::Integer, Δx::Real; derivative::Symbol = :centered)
    n ≥ 3 || throw(ArgumentError("PoissonFFT1D needs at least 3 points, got $n"))
    Δx > 0 || throw(ArgumentError("PoissonFFT1D needs Δx > 0, got $Δx"))
    derivative in (:centered, :spectral) || throw(ArgumentError(
        "derivative must be :centered or :spectral, got :$derivative"))
    h = float(Δx)
    T = typeof(h)
    k = T.(collect(2π*rfftfreq(n, 1/h)))
    spectral = derivative === :spectral
    if spectral
        coefficient = [j == 1 || (iseven(n) && j == length(k)) ? zero(Complex{T}) : -im/k[j]
                       for j in eachindex(k)]
    else
        # φ'' = ρ, the convention the centred difference below expects
        coefficient = [j == 1 ? zero(T) : -1/k[j]^2 for j in eachindex(k)]
    end
    return PoissonFFT1D{T, typeof(coefficient)}(n, h, coefficient, spectral)
end

struct PoissonWorkspace{T, P, Q}
    ρ::Vector{T}
    F::Vector{Complex{T}}
    φ::Vector{T}
    forward::P
    backward::Q
end

"""
    workspace(p::PoissonFFT1D)

The buffers and FFT plans for `p`. Planning is not thread-safe in FFTW, so
build workspaces before starting tasks, one per task.
"""
function workspace(p::PoissonFFT1D{T}, n::Integer = p.n) where {T}
    n == p.n || throw(DimensionMismatch("solver built for $(p.n) points, workspace asked for $n"))
    ρ = zeros(T, n)
    F = zeros(Complex{T}, n ÷ 2 + 1)
    φ = zeros(T, n)
    return PoissonWorkspace(ρ, F, φ, plan_rfft(ρ), plan_irfft(F, n))
end

"""
    solve!(e, ρ, p::PoissonFFT1D[, ws])

Write into `e` the field of the charge density `ρ`. Without `ws` a workspace is
allocated for the call.
"""
solve!(e, ρ, p::PoissonFFT1D) = solve!(e, ρ, p, workspace(p))

function solve!(e::AbstractVector, ρ::AbstractVector, p::PoissonFFT1D, ws::PoissonWorkspace)
    n = p.n
    length(e) == n && length(ρ) == n || throw(DimensionMismatch(
        "solver built for $n points, got e of $(length(e)) and ρ of $(length(ρ))"))
    Base.require_one_based_indexing(e, ρ)
    F, φ = ws.F, ws.φ
    # through the workspace's own array, so that any AbstractVector -- a view,
    # a column -- meets a plan made for an aligned, contiguous buffer
    copyto!(ws.ρ, ρ)
    mul!(F, ws.forward, ws.ρ)
    F[1] = zero(eltype(F))
    c = p.coefficient
    for i = 2:length(F)
        F[i] = c[i] * F[i]
    end
    if p.spectral
        mul!(φ, ws.backward, F)
        copyto!(e, φ)
    else
        mul!(φ, ws.backward, F)
        # `/(2Δx)` rather than 0.1's `0.5*…/Δx`: one rounding of the same
        # quotient, so the same bits in Float64, and no Float64 literal to
        # widen a Float32 solve.
        Δx = p.Δx
        e[1] = (φ[2] - φ[end])/(2Δx)
        for i = 2:n-1
            e[i] = (φ[i+1] - φ[i-1])/(2Δx)
        end
        e[n] = (φ[1] - φ[end-1])/(2Δx)
    end
    return e
end

"""
    generate_solver(ρ₀, Δx)

The 0.1 closure form, kept as a wrapper: `solve!(e, ρ)` over a
[`PoissonFFT1D`](@ref) with the centred derivative and one captured workspace,
so the closure is not safe to share between tasks. Deprecated; use
`PoissonFFT1D` with a `workspace` per task.
"""
function generate_solver(ρ₀, Δx)
    Base.depwarn("PoissonFourier1D.generate_solver is deprecated; use " *
                 "PoissonFFT1D(n, Δx) with workspace and solve!", :generate_solver)
    p = PoissonFFT1D(length(ρ₀), Δx)
    ws = workspace(p)
    return (e, ρ) -> solve!(e, ρ, p, ws)
end

end
