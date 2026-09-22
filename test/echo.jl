# The plasma echo: free streaming, one velocity kick, and what theory says the
# kick brings back -- exactly with the field off (`echo_closed_form`), and to
# second order with it on (`echo_second_order`).
#
# Shared by `test/VlasovSolver/test_echo.jl` and `test_verification.jl`, which
# assert it, and `verification/plasma-echo.jl`, which draws it -- for the reason
# `wakefield` is shared: the script and the tests that assert its claims run one
# setup. The self-consistent run itself, `self_consistent_echo`, is in
# `verification_harness.jl`, because it is built on `vlasov_poisson`.
#
# The files in `test/` include it as
# `@isdefined(echo_closed_form) || include(...)`, so that each still runs on its
# own while a full test run, which puts them all in `Main`, evaluates this one
# once: a repeat would redefine every method here, and `Pkg.test` runs with
# `--warn-overwrite=yes`.

using Vasilek
using SpecialFunctions: besselj1

# `Z`, `Zprime` and `dielectric`, for the second-order theory.
@isdefined(landau_root) || include(joinpath(@__DIR__, "dispersion.jl"))

"""
    echo_closed_form(t; α, ε, τ, k₁, k₂)

Complex amplitude of the `k₃ = k₂ − k₁` density mode, normalised so that
`n = a·cos k₃x` returns `a`, for free streaming from

    f₀(x, v) = (2π)^(-1/2) exp(−v²/2) · (1 + α cos k₁x)

with an instantaneous kick `v → v + ε cos k₂x` at `t = τ`:

    A₃(t) = −iα·J₁(k₃ε(t − τ))·exp(−(k₃t − k₂τ)²/2)        t > τ, zero before

that is, `n₃ = α·J₁(…)·exp(…)·sin k₃x`.

**Derivation.** Follow a particle: `x(t) = x₀ + v₀t + ε(t − τ)·cos k₂(x₀ + v₀τ)`
after the kick. Expand `exp(−ik₃x)` by Jacobi–Anger,
`exp(−iz cos θ) = Σₘ (−i)ᵐ Jₘ(z) exp(imθ)`. The average over `x₀` keeps the
terms with `mk₂ + qk₁ = k₃`, where `q ∈ {0, ±1}` runs over the Fourier content of
`f₀`; with `k₁` and `k₂` the second and third modes of the box, `3m + 2q = 1`
has the one solution `m = 1, q = −1`. The average of `exp(−iv₀(k₃t − k₂τ))` over
the Maxwellian is the Gaussian.

Nothing was linearised, so the form is exact in `α` -- free streaming is linear
in `f` -- and in `ε`, the Bessel series having been summed rather than truncated.
That difference is not small at the defaults: the linear-in-`ε` form,
`J₁(z) → z/2`, puts the peak 14.6% higher and 0.10 later.

**Where the echo is.** The Gaussian is centred where the phase-mixing exponent
vanishes, `t_e = k₂τ/(k₂ − k₁)`, the time at which the kick has unwound the shear
the seed accumulated before it. The peak of `|A₃|` comes later than `t_e`, by the
growth of `J₁` across the pulse: 15.28 against 15 at the defaults. `t_e` is set
by the physics alone and recurrence, `2π/(kΔv)`, by the grid alone; the defaults
keep the second (62.8 for the seeded mode at `Δv = 0.1`, 31.4 at 0.2) outside
the window.
"""
function echo_closed_form(t; α, ε, τ, k₁, k₂)
    t ≤ τ && return zero(ComplexF64)
    k₃ = k₂ - k₁
    return -im*α*besselj1(k₃*ε*(t - τ))*exp(-(k₃*t - k₂*τ)^2/2)
end

"""
    ballistic_echo(scheme_x, scheme_v; L = 4π, m₁ = 2, m₂ = 3, Nx = 64, Δv = 0.1,
                   vmax = 6.0, Δt = 0.02, α = 0.1, ε = 0.2, τ = 5.0, tmax = 22.0,
                   kick_courant = 0.5, snapshots = ())

Run the setup of [`echo_closed_form`](@ref): free streaming in `x`, one row per
velocity as in `test_free_streaming.jl`, and at `t = τ` the kick, an advection in
`v` of every column by `ε cos k₂x`. `k₁` and `k₂` are the `m₁`-th and `m₂`-th
modes of the box. Either scheme may be given as a function of the initial `f`
returning one -- `f -> PFC(fmin = 0.0, fmax = maximum(f))` -- which is how `PFC`
is bounded by the distribution it will carry.

Returns `t`; `k = (k₁, k₂, k₃)`; `modes`, the `Nt × 3` matrix of complex density
amplitudes at those three wavenumbers; `kick`, the index of `τ` in `t`; the grids,
the final `f`, and copies of `f` at the `snapshots` times.

**Row `s` of `modes` is taken at `t[s]`, before the step, and the row at `τ`
before the kick** -- which makes that sample the seeded mode at its most hidden,
and a snapshot at `τ` the distribution the kick is handed.

**The kick is applied in sub-steps** of Courant number at most `kick_courant`,
because every scheme here except the semi-Lagrangian one is restricted to
`|c| ≤ 1`. At fixed `x` successive shifts compose exactly, so the sub-steps are
the same instantaneous kick. What they change is how many times a scheme
interpolates the filaments, which is part of what is being measured, and so the
count is set by the grid rather than by the scheme: four sub-steps of Courant up
to 0.5 at `ε = 0.2`, `Δv = 0.1`, for every scheme.
"""
function ballistic_echo(scheme_x, scheme_v; L = 4π, m₁ = 2, m₂ = 3, Nx = 64, Δv = 0.1,
                        vmax = 6.0, Δt = 0.02, α = 0.1, ε = 0.2, τ = 5.0, tmax = 22.0,
                        kick_courant = 0.5, snapshots = ())
    k₁, k₂ = 2π*m₁/L, 2π*m₂/L
    k = (k₁, k₂, k₂ - k₁)
    Δx = L/Nx
    x = [(j-1)*Δx for j = 1:Nx]
    v = collect(-vmax:Δv:vmax)
    Nv = length(v)
    f = [exp(-v[i]^2/2)/sqrt(2π)*(1 + α*cos(k₁*x[j])) for i = 1:Nv, j = 1:Nx]

    t = collect(0:Δt:tmax)
    # The kick and every snapshot have to land on a step of the run. A snapshot
    # that does not is never written, and `similar` would hand back whatever the
    # allocator left there -- measured 1.1e293 for a snapshot past `tmax` -- as
    # if it were a distribution function.
    on_step(ts, what) = begin
        s = round(Int, ts/Δt) + 1
        1 ≤ s ≤ length(t) && isapprox(t[s], ts; atol = 1e-9*Δt) ||
            error("$what at t = $ts does not fall on a step of Δt = $Δt in [0, $(t[end])]")
        s
    end
    kick = on_step(τ, "the kick")
    nsub = max(1, ceil(Int, abs(ε)/(kick_courant*Δv)))
    shots = [similar(f) for _ in snapshots]
    shot_at = [on_step(ts, "a snapshot") for ts in snapshots]

    modes = zeros(ComplexF64, length(t), 3)
    scheme_x = scheme_x isa AbstractAdvection1D ? scheme_x : scheme_x(f)
    scheme_v = scheme_v isa AbstractAdvection1D ? scheme_v : scheme_v(f)
    wsx, wsv = workspace(scheme_x, Nx), workspace(scheme_v, Nv)
    bufx, bufv, col = zeros(Nx), zeros(Nv), zeros(Nv)
    n = zeros(Nx)
    for s in eachindex(t)
        for j = 1:Nx
            n[j] = Δv*sum(view(f, :, j))
        end
        for (m, km) in enumerate(k)
            modes[s, m] = 2*sum(n[j]*cis(-km*x[j]) for j = 1:Nx)/Nx
        end
        for (i, si) in enumerate(shot_at)
            si == s && copyto!(shots[i], f)
        end
        s == length(t) && break

        if s == kick
            for j = 1:Nx
                c = ε*cos(k₂*x[j])/(nsub*Δv)
                copyto!(col, view(f, :, j))
                for _ = 1:nsub
                    advect!(bufv, col, scheme_v, c, wsv)
                    copyto!(col, bufv)
                end
                f[:, j] .= col
            end
        end
        for i = 1:Nv
            advect!(bufx, view(f, i, :), scheme_x, v[i]*Δt/Δx, wsx)
            f[i, :] .= bufx
        end
    end
    return (; t, k, modes, kick, x, v, f, nsub, snapshots = shots)
end

# ------------------------------------------------ with the plasma's own field

"""
    resonant_dielectric(k, v)

`D(k, ω)` of a unit Maxwellian and `∂D/∂ω`, both on the real axis at `ω = kv`, the
frequency a particle of velocity `v` sees a wave of wavenumber `k` at:

    D = 1 + (1 + ζZ(ζ))/k²,    ∂D/∂ω = (Z(ζ) + ζZ′(ζ))/(√2·k³),    ζ = v/√2

`D` is `dielectric(k*v, k, maxwellian)` -- `test_echo.jl` checks that -- written
out because the derivative is needed beside it.
"""
function resonant_dielectric(k, v)
    ζ = v/sqrt(2)
    return 1 + (1 + ζ*Z(ζ))/k^2, (Z(ζ) + ζ*Zprime(ζ))/(sqrt(2)*k^3)
end

"""
    echo_filament(v; α, k₁)

The ballistic part of a seeded mode once its own field has Landau-damped, and its
`v`-derivative: the seed leaves `f_{k₁} → h₁(v)·exp(−ik₁vt)`, with

    h₁(v) = (α/2)·M(v)·(1 + 1/k₁²)/D(k₁, k₁v)

The general form is `(α/2)M − M′·ê₁(k₁v)`, `ê₁` being the Laplace transform of
the seed's field at the particle's resonance. For a Maxwellian, `M′ = −vM`, and
the identity `D = 1 + (1 − ωW)/k²` for `W = ∫M/(ω − kv)dv` collapses it to the
line above. `test_echo.jl` checks it against the general form with `ê₁` taken
from a time-domain solve of the seed problem, which shares nothing with this
route but the Maxwellian.

Without the field the filament would be `(α/2)M` exactly, and the difference is
not a correction. At `k₁ = 1` the two agree at `v = 0`, where `D = 1 + 1/k₁²`
exactly; `|h₁|` is 1.35 times `(α/2)M` at `v = 1`, 2.51 at `v = 2`, and 2.66 at
`v = 2.34`, beside the phase velocity 2.05 of the mode the seed launched. That is
where the Landau-damped field put its energy.
"""
function echo_filament(v; α, k₁)
    D, Dω = resonant_dielectric(k₁, v)
    c = (α/2)*(1 + 1/k₁^2)*exp(-v^2/2)/sqrt(2π)
    return c/D, c*(-v/D - k₁*Dω/D^2)
end

"""
    maxwellian_response(S, t, k)

Total density at wavenumber `k` of a unit-Maxwellian plasma in which the density
`S(t)` streams ballistically, on the uniform range `t`: the solution of

    n(t) = S(t) − ∫ (t − t′)·exp(−k²(t − t′)²/2)·n(t′) dt′      over t′ ∈ [t₁, t]

The kernel is what a field `e = n/(ik)` does to the Maxwellian,
`−∫M′(v)·exp(−ikv(t − t′))dv = −ik(t − t′)·exp(−k²(t − t′)²/2)`, so this is
`n = S/D(k, ω)` in the time domain. For a cold plasma the kernel is `t − t′`
and the equation is `n″ = S″ − n`: an oscillation at the plasma frequency, which
fixes the sign. Trapezoidal, and the diagonal term vanishes with the kernel, so
it is explicit.
"""
function maxwellian_response(S, t::AbstractRange, k)
    Δt = step(t)
    K = [(m*Δt)*exp(-(k*m*Δt)^2/2) for m in 0:length(t)-1]
    n = similar(S)
    for i in eachindex(t)
        acc = zero(eltype(S))
        for j in 1:i-1
            acc += (j == 1 ? 0.5 : 1.0)*K[i-j+1]*n[j]
        end
        n[i] = S[i] - Δt*acc
    end
    return n
end

"""
    echo_second_order(t; α, ε, τ, k₁, k₂, v = -9.0:0.002:9.0, field = true)

The echo of [`echo_closed_form`](@ref) with the plasma's own field on: second
order -- first in `α`, first in `ε` -- and long after each perturbation compared
with the damping of the field it launched. Returns `(; echo, source)`, the density
amplitude of the `k₃ = k₂ − k₁` mode and the ballistic density it answers, both
normalised like `echo_closed_form`, on `t`, a uniform range that starts after `τ`.

Three linear responses of the Maxwellian, composed; `D` is
[`resonant_dielectric`](@ref).

1. **The seed's filament**, [`echo_filament`](@ref): `h₁(v)·exp(−ik₁vt)`.
2. **The kick, meeting it.** The shift `v → v + ε cos k₂x` turns the `−k₁` half of
   the filament into the `k₃` mode, and so does the field the kick itself induces
   at `k₂`, for as long as it lasts. After both,

       f₃ ⊃ −Q(v)·exp(−iv(k₃t − k₂τ)),
       Q = (ε/2)·[(ik₁τ·h₁* + h₁*′)/D(k₂, k₂v) − k₁·h₁*·∂D/∂ω(k₂, k₂v)/D(k₂, k₂v)²]

   The first term is the kick, screened: the induced field adds
   `(ε/2)(1/D − 1)` to the external `ε/2`. The second is the induced field's
   *duration*, during which the filament keeps shearing -- the first moment in
   time of that field, which is what the `ω`-derivative of its transform is.
3. **The echo's own field.** The ballistic density
   `S(t) = −∫Q·exp(−iv(k₃t − k₂τ))dv` polarises the plasma at `k₃`, and the
   total is [`maxwellian_response`](@ref) to it.

`field = false` drops all three -- `h₁ = (α/2)M`, `D = 1`, no response -- and must
give the small-`ε` limit of `echo_closed_form`, which is how the signs were
checked; `test_echo.jl` asserts that too.

**What is left out.** What remains of the seed's own field when the kick arrives,
`exp(−γ(k₁)τ)` of it -- 2.0e-4 at `τ = 10` against 1.4e-2 at `τ = 5`, which is
why the self-consistent run waits -- and everything that field does afterwards;
terms of order `α²`, which moved the comparison with the run by 3e-5 of the peak
when `α` was cut from 0.01 to 0.001; and the curvature of `J₁`, 1.2e-3 at
`ε·k₁τ = 0.1`, which is why `ε` is small. The residual against the run is then
the run's own truncation, and `test_verification.jl` shows it converging.
"""
function echo_second_order(t::AbstractRange; α, ε, τ, k₁, k₂, v = -9.0:0.002:9.0,
                           field = true)
    first(t) > τ || error("the echo comes after the kick; t starts at $(first(t)) ≤ τ = $τ")
    k₃ = k₂ - k₁
    Q = map(v) do u
        if field
            h, h′ = echo_filament(u; α, k₁)
            D, Dω = resonant_dielectric(k₂, u)
            (ε/2)*((im*k₁*τ*conj(h) + conj(h′))/D - k₁*conj(h)*Dω/D^2)
        else
            (ε/2)*(α/2)*(im*k₁*τ - u)*exp(-u^2/2)/sqrt(2π)
        end
    end
    S = [-step(v)*sum(Q[j]*cis(-v[j]*(k₃*s - k₂*τ)) for j in eachindex(v)) for s in t]
    n = field ? maxwellian_response(S, t, k₃) : S
    return (echo = 2 .* n, source = 2 .* S)
end
