# Shared 1D1V Vlasov–Poisson driver for the extended verification tests.
# Deliberately close to what the verification notebooks do, so the two agree.

using Vasilek
using Vasilek: StrangSplitting, FDTD1D, PoissonFourier1D
using NumericalIntegration, FFTW

# The kinetic dispersion relation: `landau_root`, `two_stream_warm` and the
# pieces they are built from. Separate file because nothing in it runs a
# simulation, and `test_dispersion.jl` exercises it without the Strang loop.
include(joinpath(@__DIR__, "dispersion.jl"))

# The plasma echo's closed form, its second-order theory, and the free-streaming
# run; `self_consistent_echo` below is the same experiment through
# `vlasov_poisson`.
include(joinpath(@__DIR__, "echo.jl"))

"Spectral Poisson solve on a uniform x grid, e = -dφ/dx with φ'' = -ρ."
function make_poisson(x)
    Δx = x[2] - x[1]
    L = x[end] - x[1] + Δx
    ω = 2π*collect(0.0:1.0/L:0.5/Δx)
    ω[1] = ω[2]
    return function (e, ρ)
        F = FFTW.rfft(ρ)
        φ = FFTW.irfft(F ./ (-ω.^2), length(ρ))
        e[:] = vcat(0.5*[φ[2] - φ[end]],
                    0.5*(φ[3:end] - φ[1:end-2]),
                    0.5*[φ[1] - φ[end-1]]) ./ Δx
        return e
    end
end

cell_widths(z) = vcat([z[2]-z[1]], 0.5*(z[3:end] - z[1:end-2]), [z[end]-z[end-1]])

"""
    nlogn(u)

`-u·log u`, and zero at `u ≤ 0`, for the entropy integrand.

A scalar function rather than an `ifelse` inside the broadcast, because
`ifelse` is an ordinary call and evaluates **both** arguments: written that way
the guard does not guard, and `log` is handed the negative value anyway. That is
not hypothetical here -- `LaxWendroff` and cubic `SemiLagrangian` drive `f` to
-0.094 and -0.098 on the large-amplitude case in
`verification/scheme-comparison.jl`, against a peak of 0.6, and the entropy
diagnostic threw `DomainError` on both until this was split out. That is a 16%
undershoot of the peak rather than round-off leaking below zero, which is worth
stating precisely: it is the size of the overshoot that makes the guard a
statement about the schemes rather than about floating point.
"""
nlogn(u) = u > 0 ? -u*log(u) : zero(u)

"""
    line_advector(scheme, Δz)

Wrap `scheme` as an in-place `(column, α)` advector, where `α` is always a
**displacement** -- a length -- whatever the scheme's own fourth argument means.

`PFCNonUniform` takes a displacement already. Every other scheme takes a Courant
number, which only exists on a uniform grid, so the wrapper divides by the
spacing for those and **refuses** a non-uniform grid rather than picking one of
its spacings and being quietly wrong by the ratio between them. That asymmetry
is a documented wart of the advection API (see `docs/normalization.md`); this is
the one place the verification runs have to absorb it.
"""
function line_advector(scheme::PFCNonUniform, Δz)
    n = length(Δz)
    ws = workspace(scheme, n)
    buf = Vector{Float64}(undef, n)
    return (col, α) -> (advect!(buf, col, scheme, α, ws); copyto!(col, buf))
end

function line_advector(scheme, Δz)
    n = length(Δz)
    h = Δz[1]
    all(d -> isapprox(d, h; rtol = 1e-12), Δz) || error(
        "$(nameof(typeof(scheme))) takes a Courant number, which a non-uniform grid " *
        "does not have (spacings range over $(extrema(Δz))); use PFCNonUniform here")
    ws = workspace(scheme, n)
    buf = Vector{Float64}(undef, n)
    return (col, α) -> (advect!(buf, col, scheme, α/h, ws); copyto!(col, buf))
end

"""
    vlasov_poisson(x, v, f₀, t; scheme_x, scheme_v, invariants = false, modes = (),
                   nᵢ = nothing, renormalize = nᵢ === nothing, stop_at_courant = false)

Strang-split Vlasov–Poisson on a static grid.

Returns a NamedTuple with `ε_e` (electric energy) and `ε` (total energy)
histories. With `invariants = true` it also returns the `mass`, `momentum`,
`l2` and `entropy` histories -- the conserved quantities that are *not* the
energy, and that nothing asserted until now. `modes = (k₁, k₂, …)` adds
`E_modes`, an `Nt × length(modes)` matrix of complex field amplitudes from
[`mode_amplitude`](@ref), which is what separates a mode from its harmonics
where `ε_e` cannot.

**Row `k` of `E_modes` (and of `ε_e`) is sampled at `t[k] + Δt/2`, not `t[k]`.**
The field is solved inside the step, after Strang's first `x` half-step, and
recorded from there. A rate or a frequency cannot see a constant time offset;
a *phase* can, and does: removing a Doppler factor `exp(-ikut)` at `t[k]` leaves
`k·u·Δt/2` behind on every sample, which the Galilean test in
`test_verification.jl` once reported as the grid's own non-invariance.

**The invariants use the cell-width sum `Σ f ΔvΔx`, not `integrate`.** That is
the quadrature the schemes actually conserve: `PFC` is a flux form, so what
leaves one cell enters its neighbour and the full-weight sum is preserved
exactly. The trapezoid halves the two endpoint weights, which no flux
conservation law protects, and measuring with it reports a drift that belongs to
the quadrature rather than to the scheme. Measured over the k = 0.5 Landau run,
875 steps: mass drifts 2.8e-16 by the cell-width sum against **2.4e-4** by the
trapezoid, and momentum stays at 5.3e-16 against 7.3e-4. Both trapezoid figures
are the endpoint weighting, not the solver.

The energy histories above keep `integrate`, because they are compared with
tolerances of half a percent where the difference is irrelevant, and because
changing them would silently move numbers the notebooks quote.

**The ions are a fixed background**, a Maxwellian's density on the grid unless
`nᵢ` gives a profile over `x`. Without `nᵢ`, `f` is rescaled on entry so that the
two integrate to the same charge; with it, `f` is taken as it comes; and
`renormalize` overrides either. The rescaling is the trapezoid over `x` of both,
which on a periodic grid weights the two end points by half, so it is exact only
while `nₑ` and `nᵢ` are proportional -- a uniform background, which is what a
caller that does not pass `nᵢ` gets. One that does is handing over a matched
pair, and the rescaling would unmatch it: for the equilibrium of
[`bgk_equilibrium`](@ref) it is 1 − 1.9e-3, and applying it doubles the
equilibrium's departure from itself through `t = 50`, from 3.85e-3 to 8.07e-3 in
the field and from 1.52e-3 to 3.17e-3 in `f`. That is inside the tolerances the
equilibrium is held to, so nothing downstream would catch a caller who forgot to
switch it off; passing `nᵢ` switches it off instead.

`scheme_x` and `scheme_v` default to `PFCNonUniform` on the two grids, which is
what the verification notebooks use and what every previous caller got. They are
arguments so that the same driver can measure what the physics costs under a
*different* scheme, which is what `verification/scheme-comparison.jl` does, and
so that a refinement study can hold the scheme fixed while moving the grid.

**Either may also be a function of the initial `f` that returns a scheme**, as
in `f -> PFC(fmin = 0.0, fmax = maximum(f))`, and is called with the `f` the run
actually starts from. That is the only way for a caller to bound a scheme by
that distribution: the driver rescales what it is handed before it runs -- by
0.8% at α = 0.5 -- so a bound taken from `f₀` beforehand sits below the rescaled
maximum and trips `PFC`'s own check on the first call.

**The defaults' bounds are those of `f` on entry: 0 and its maximum.** By
Liouville's theorem the exact solution keeps both, and PFC's limiter exists to
keep a run between the bounds it is given, so they belong to the run rather than
to the driver. They were the constant 1, which nothing chose and `PFCNonUniform`
did not then check. Above it the limiter's `2(fmax − f)` goes negative and the
scheme corrupted the run without a word: an equilibrium whose trapped population
peaks at 1.79 ended 44% of its peak away from itself with the bound at 1 -- a run
the scheme now refuses on its first call. Below it, where every run until then
sat, the bound never engaged at all.

Tight, it does engage, at the maximum, and that has a measured price: the limiter
clips the reconstruction in the peak cell, and the α = 0.05 round trip in
`test_verification.jl` converges at second order (×4.4, then ×4.1) where it
converged at third (×5.7, ×7.1). Elsewhere the numbers moved little and are
updated where they are quoted. The `fmax` history of a Landau, a strong Landau,
a two-stream and an equilibrium run tops out exactly at the bound, or below, to
the last bit. One run did leave it, and recorded no history to show it: the
`a = 0.6` two-stream case, past its velocity Courant limit. `PFCNonUniform` now
checks every call, so none can leave it quietly, and that run stops at the limit
instead.

**`stop_at_courant = true` ends the run at the velocity Courant limit.** The
x-sweeps move `f` by `vΔt/2`, which the grid fixes; the v-sweep moves it by
`eΔt`, which grows with the field, and past one cell per step PFC's flux reaches
beyond the neighbouring cell and the scheme leaves the bounds it is built on --
which, checked, it refuses a few steps later. With the flag the loop ends after
the first step whose `max|e|Δt/min Δv` exceeds 1, and every history ends with
that step: `length(ε_e) - 1` steps ran, and `f` is the state after the last.
The unstable runs need it; see [`two_stream`](@ref).
"""
function vlasov_poisson(x, v, f₀, t;
                        scheme_x = nothing, scheme_v = nothing, invariants = false,
                        modes = (), nᵢ = nothing, renormalize = nᵢ === nothing,
                        stop_at_courant = false)
    Δx = cell_widths(x)
    Δv = cell_widths(v)

    if nᵢ === nothing
        fᵢ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. Δx/Δx)'
        nᵢ = integrate(v, fᵢ)
    else
        length(nᵢ) == length(x) || throw(DimensionMismatch(
            "nᵢ has $(length(nᵢ)) points, the x grid $(length(x))"))
        nᵢ = collect(float.(nᵢ))
    end
    Nᵢ = integrate(x, nᵢ)

    f = copy(f₀)
    renormalize && (f .*= Nᵢ/integrate(x, integrate(v, f)))
    g = f'

    bound = maximum(f)
    pick(s, widths) = s === nothing ? PFCNonUniform(widths; fmin = 0.0, fmax = bound) :
                      s isa AbstractAdvection1D ? s : s(f)
    sx, sv = pick(scheme_x, Δx), pick(scheme_v, Δv)

    advect_x! = line_advector(sx, Δx)
    advect_v! = line_advector(sv, Δv)

    solve_poisson! = make_poisson(x)
    e = similar(x)
    ε_e = similar(t)
    ε = similar(t)
    # One scratch matrix, reused: the invariants are five integrals of the same
    # shape as `f`, and allocating a temporary per integral per step dominated
    # the step itself when this was first written.
    tmp = similar(f)
    wt       = invariants ? Δv .* Δx' : nothing      # the flux-form quadrature
    mass     = invariants ? similar(t) : nothing
    momentum = invariants ? similar(t) : nothing
    l2       = invariants ? similar(t) : nothing
    entropy  = invariants ? similar(t) : nothing
    # `f` is a distribution function: negative values are unphysical, and a
    # scheme that produces them can be accurate on a linear diagnostic while
    # being unusable on a nonlinear run. Tracked because the comparison study
    # would otherwise rank a non-positive scheme first without saying so.
    fmin     = invariants ? similar(t) : nothing
    # And the maximum, which Liouville's theorem forbids to grow just as it
    # forbids the minimum to fall; the defaults hold their limiter to the
    # initial one, and this is how a run shows whether it stayed there.
    fmax     = invariants ? similar(t) : nothing
    # Per-mode field amplitudes, when asked for. `ε_e` sums every mode in the
    # box, which is fine while one of them dominates and misleading the moment
    # another does -- the recurrence of the second harmonic arrives at half the
    # time the seeded mode's does, and in `ε_e` it is indistinguishable from the
    # seeded mode coming back early.
    E_modes  = isempty(modes) ? nothing : zeros(ComplexF64, length(t), length(modes))

    steps = length(t) - 1
    for k in 1:length(t)-1
        Δt = t[k+1] - t[k]
        vΔt(_) = v*Δt
        function eΔt(ff)
            nₖ = vec(sum(ff'.*Δv, dims = 1))
            solve_poisson!(e, nₖ - nᵢ)
            return e*Δt
        end
        StrangSplitting.make_time_step_2d!((g, f), (vΔt, eΔt), (advect_x!, advect_v!))
        if E_modes !== nothing
            for (j, km) in enumerate(modes)
                E_modes[k, j] = mode_amplitude(e, x, km)
            end
        end
        ε_e[k] = integrate(x, e.^2)
        @. tmp = f*v^2
        ε[k] = integrate(x, integrate(v, tmp)) + ε_e[k]

        if invariants
            @. tmp = f*wt
            mass[k] = sum(tmp)
            @. tmp = f*v*wt
            momentum[k] = sum(tmp)
            @. tmp = f*f*wt
            l2[k] = sum(tmp)
            @. tmp = nlogn(f)*wt
            entropy[k] = sum(tmp)
            fmin[k] = minimum(f)
            fmax[k] = maximum(f)
        end

        # This step's v-sweep moved `f` by `e*Δt`; the one that crossed a cell
        # is the last.
        if stop_at_courant && maximum(abs, e)*Δt > minimum(Δv)
            steps = k
            break
        end
    end
    for h in (ε_e, ε, mass, momentum, l2, entropy, fmin, fmax)
        h === nothing || (h[steps+1] = h[steps])
    end
    E_modes === nothing || (E_modes[steps+1, :] = E_modes[steps, :])
    if steps < length(t) - 1
        ε_e, ε, mass, momentum, l2, entropy, fmin, fmax =
            map(h -> h === nothing ? nothing : h[1:steps+1],
                (ε_e, ε, mass, momentum, l2, entropy, fmin, fmax))
        E_modes = E_modes === nothing ? nothing : E_modes[1:steps+1, :]
    end
    # `f` comes back too. It costs nothing -- the array exists either way -- and
    # it is the only way to ask a question about the distribution rather than
    # about a moment of it, which is what the reversibility test needs.
    return (; ε_e, ε, mass, momentum, l2, entropy, fmin, fmax, E_modes, f)
end

"""
    self_consistent_echo(; L = 4π, m₁ = 2, m₂ = 3, Nx = 128, Δv = 0.05, vmax = 6.0,
                         Δt = 0.02, α = 0.01, ε = 0.01, τ = 10.0, tmax = 40.0,
                         kick_courant = 0.5)

[`ballistic_echo`](@ref) with the field on: `vlasov_poisson` from `t = 0` to `τ`,
the kick, and `vlasov_poisson` again from `τ` on the kicked distribution -- a
restart the driver supports as it stands, since it takes whatever `f` it is
handed. The kick is `ballistic_echo`'s, sub-steps included, through the driver's
own velocity scheme.

Returns `t`, `k = (k₁, k₂, k₃)` and `modes`, the density amplitudes `ik·E_k` at
those wavenumbers -- densities rather than fields, so that they compare with
`echo_second_order` and `echo_closed_form` directly. **`t` holds the times the
modes were sampled at, which are mid-step, `t + Δt/2`** (see `vlasov_poisson`);
the row the driver duplicates at the end of each segment is dropped.

**The defaults are not `ballistic_echo`'s**, and each difference is there for the
theory. `echo_second_order` is second order, so both amplitudes are small,
`α = ε = 0.01`. It holds once the seed's own field has damped, and at `k₁ = 1`
that field damps at `γ = 0.851`: waiting until `τ = 10` leaves 2.0e-4 of it where
`τ = 5` would leave 1.4e-2. The filament the kick then meets is twice as fine,
`2π/(k₁τ) = 0.63`, which `Δv = 0.05` resolves with thirteen cells, and the echo
arrives at `t_e = 30`. And `Nx = 128`, because the x-sweep is what limits the
comparison: the run is 1.6e-2 of the peak from the theory at `Nx = 64` and 4.8e-3
at 128, where halving `Δv` or `Δt` at `Nx = 64` instead leaves it at 1.5e-2 and
1.7e-2.
"""
function self_consistent_echo(; L = 4π, m₁ = 2, m₂ = 3, Nx = 128, Δv = 0.05, vmax = 6.0,
                              Δt = 0.02, α = 0.01, ε = 0.01, τ = 10.0, tmax = 40.0,
                              kick_courant = 0.5)
    k₁, k₂ = 2π*m₁/L, 2π*m₂/L
    k = (k₁, k₂, k₂ - k₁)
    Δx = L/Nx
    x = [(j-1)*Δx for j = 1:Nx]
    v = collect(-vmax:Δv:vmax)
    f₀ = [exp(-u^2/2)/sqrt(2π)*(1 + α*cos(k₁*y)) for u in v, y in x]

    before = collect(0:Δt:τ)
    isapprox(before[end], τ; atol = 1e-9*Δt) ||
        error("the kick at τ = $τ does not fall on a step of Δt = $Δt")
    seed = vlasov_poisson(x, v, f₀, before; modes = k)

    f = copy(seed.f)
    widths = cell_widths(v)
    kick! = line_advector(PFCNonUniform(widths; fmin = 0.0, fmax = maximum(f)), widths)
    nsub = max(1, ceil(Int, abs(ε)/(kick_courant*Δv)))
    for j in eachindex(x), _ = 1:nsub
        kick!(view(f, :, j), ε*cos(k₂*x[j])/nsub)
    end

    after = collect(τ:Δt:tmax)
    echo = vlasov_poisson(x, v, f, after; modes = k)

    t = vcat(before[1:end-1], after[1:end-1]) .+ Δt/2
    E = vcat(seed.E_modes[1:end-1, :], echo.E_modes[1:end-1, :])
    return (; t, k, modes = E .* transpose(im .* collect(k)), x, v)
end

# --------------------------------------------------- nonlinear equilibria

"""
    bgk_distribution(ψ; T_trapped = 1.0)

`F(W)` of a stationary solution in the potential energy `U = −ψ cos kx`, as a
function of the particle energy `W = v²/2 + U`: Maxwellian for the passing
particles, `W ≥ ψ`, and at temperature `T_trapped` below the separatrix, the two
joined continuously at `W = ψ`.

`T_trapped = 1` is the Maxwell–Boltzmann equilibrium `M(v)·exp(ψ cos kx)`,
analytic across the separatrix. Any other value is a BGK mode proper -- a trapped
population the passing one does not determine -- and puts a kink in `F` exactly
on the separatrix, which is the structure it exists to test.
"""
bgk_distribution(ψ; T_trapped = 1.0) =
    W -> W ≥ ψ ? exp(-W)/sqrt(2π) : exp(-ψ - (W - ψ)/T_trapped)/sqrt(2π)

"""
    bgk_equilibrium(; ψ = 0.5, k = 0.5, T_trapped = 1.0, Nx = 64, Δv = 0.1,
                    vmax = 6.0, Δt = 0.05, tmax = 50.0, ψ_ions = ψ, ion_sign = +1)

Run a nonlinear equilibrium: `f₀ = F(v²/2 − ψ cos kx)` with `F` from
[`bgk_distribution`](@ref), over one wavelength, on the ion background that
holds it still.

Any function of the energy is a stationary solution of the Vlasov equation; what
makes it one of Vlasov--Poisson is a charge density whose field is the force
`−∂ₓU = −ψk sin kx`. With the uniform ions every other run uses there is none,
so the ions are built for it: `∂ₓE = nₑ − nᵢ` in this code's convention gives

    nᵢ = nₑ + ∂ₓ²U = nₑ + ψk² cos kx

with `nₑ` the code's own moment of `f₀`. `ψ_ions` and `ion_sign` build them for a
different `ψ`, or with the sign of the Poisson equation reversed -- the two ways
of handing the run a state that is not an equilibrium.

Returns the grids, `f₀`, the ions `nᵢ` and the final `f`, the field's `k` mode
history `E` (mid-step, see `vlasov_poisson`) against `E₀`, the equilibrium's
own: `iψk` times `sin(kΔx)/(kΔx)`, the centred difference `docs/normalization.md`
documents. Also the separatrix `v_sep(x) = √(2ψ(1 + cos kx))` and the `l2` and
`entropy` histories.

**The defaults.** `ψ = 0.5` traps everything below `|v| = 1.41` at the bottom of
the well -- 85% of the particles there, two thirds of all of them -- and the
deepest bounce at `ω_B = k√ψ = 0.354`, so `t = 50` is 2.8 of their periods. The peak of `f₀` is
0.66 for `T_trapped = 1`, 0.40 for `T_trapped = 2` and 1.79 for `T_trapped = 1/2`.
"""
function bgk_equilibrium(; ψ = 0.5, k = 0.5, T_trapped = 1.0, Nx = 64, Δv = 0.1,
                         vmax = 6.0, Δt = 0.05, tmax = 50.0, ψ_ions = ψ, ion_sign = +1)
    Δx = 2π/k/Nx
    x = [(j-1)*Δx for j = 1:Nx]
    v = collect(-vmax:Δv:vmax)
    moment(F, ψ) = [Δv*sum(F(u^2/2 - ψ*cos(k*y)) for u in v) for y in x]
    F = bgk_distribution(ψ; T_trapped)
    f₀ = [F(u^2/2 - ψ*cos(k*y)) for u in v, y in x]
    nᵢ = moment(bgk_distribution(ψ_ions; T_trapped), ψ_ions) .+
         ion_sign*ψ_ions*k^2 .* cos.(k .* x)

    t = collect(0:Δt:tmax)
    r = vlasov_poisson(x, v, f₀, t; nᵢ, modes = (k,), invariants = true)
    return (; x, v, t = t[1:end-1] .+ Δt/2, f₀, nᵢ, f = r.f, E = r.E_modes[1:end-1, 1],
            E₀ = im*ψ*k*sin(k*Δx)/(k*Δx), v_sep = @.(sqrt(2ψ*(1 + cos(k*x)))),
            l2 = r.l2[1:end-1], entropy = r.entropy[1:end-1])
end

# --------------------------------------------------------- mode fitting
#
# The electric energy of a single Landau mode goes as
#
#     ε_e(t) ∝ exp(-2γt)·cos²(ω t + φ)
#
# so both the damping rate and the real frequency are recoverable from it --
# but only if the `cos²` is handled rather than ignored.

"""
    mode_amplitude(e, x, k)

Complex amplitude of the `exp(ikx)` component of `e` on the uniform grid `x`,
normalised so that `e = A·cos(kx)` returns `A`.

`ε_e` is the only field diagnostic the driver kept until now, and it is a sum
over every mode in the box. That is enough while one mode dominates and
actively misleading when a second one arrives: the recurrence of a nonlinearly
generated harmonic lands at `2π/(mkΔv)`, i.e. earlier than the seeded mode's own
by the harmonic number, and in `ε_e` the two are the same bump.
"""
mode_amplitude(e, x, k) = 2*sum(e[j]*cis(-k*x[j]) for j in eachindex(x))/length(x)

"""
    local_extrema(t, y; tmin, tmax, maxima)

Indices of the strict interior local maxima (`maxima = true`) or minima of `y`,
restricted to `tmin ≤ t[i] ≤ tmax`.
"""
function local_extrema(t, y; tmin, tmax, maxima::Bool)
    return [i for i in 2:length(y)-1 if tmin ≤ t[i] ≤ tmax &&
            (maxima ? (y[i] > y[i-1] && y[i] > y[i+1])
                    : (y[i] < y[i-1] && y[i] < y[i+1]))]
end

"""
    damping_rate(t, ε_e; tmin, tmax)

Landau damping rate `γ`, from a least-squares fit of `log ε_e` against `t`
through the **local maxima only**.

Fitting every sample instead -- which is what this suite did until now -- fits
`log(exp(-2γt)·cos²(ωt+φ))`, and `log cos²` has a pole at every null of the
oscillation. The result is dominated by how close the window edges happen to
land to a null, and it moves discontinuously as the window is nudged.

Measured at k = 0.5, where the tabulated root is 0.15336. Fitting every sample:
0.14837 to 0.15755 as the window varies over plausible choices, a spread of 6.2%
of the value -- and the window this file used, `t ∈ [5.9, 29.9]`, began *exactly*
on a minimum, which is the entire reason it reported 0.1498 (2.3% low). Moving
the start one step, to 6.0, gives 0.1532 (0.1%) from the same data.

Through the maxima the same sweep gives 0.15451 to 0.15571, a spread of 0.8%,
consistently about 1% above the analytic value. That residue is numerical
damping and does not move under refinement; the 6.2% was an artefact of the
estimator.
"""
function damping_rate(t, ε_e; tmin, tmax)
    p = local_extrema(t, ε_e; tmin = tmin, tmax = tmax, maxima = true)
    length(p) ≥ 3 || error("damping_rate needs at least 3 maxima in [$tmin, $tmax], found $(length(p))")
    A = hcat(ones(length(p)), t[p])
    return -(A \ log.(ε_e[p]))[2]/2, length(p)
end

"""
    oscillation_frequency(t, ε_e; tmin, tmax)

Real frequency `ω`, from the mean spacing of the minima of `ε_e`.

The minima are spaced by `π/ω` rather than `2π/ω`: the energy carries `cos²`,
which has twice the frequency of the field. They are used in preference to the
maxima because a null is a sharp feature whose location is well defined even
once the amplitude has decayed by several orders of magnitude, and only the
first and last are needed, so the estimate improves with the length of the
window rather than degrading with it.
"""
function oscillation_frequency(t, ε_e; tmin, tmax)
    m = local_extrema(t, ε_e; tmin = tmin, tmax = tmax, maxima = false)
    length(m) ≥ 2 || error("oscillation_frequency needs at least 2 minima in [$tmin, $tmax], found $(length(m))")
    return π/((t[m[end]] - t[m[1]])/(length(m) - 1)), length(m)
end

# ----------------------------------------------------- laser wakefield
#
# The electromagnetic study, extracted from `verification/wakefield.jl` so that
# the script and the test that asserts its claims run the same code. The script
# now calls this and plots what comes back; before, there was nothing to call.

"""
    ddx!(d, u, Δx)

`∂u/∂x` by the same centred difference, periodic in `x`, that
`PoissonFourier1D` takes to get `e` from `φ`.

The stencil is not an arbitrary choice of second-order rule: the ponderomotive
force and the electrostatic one are added together before the momentum push, so
differencing them alike is what keeps the sum from carrying a gradient that
belongs to neither.
"""
function ddx!(d, u, Δx)
    d[1] = 0.5*(u[2] - u[end])/Δx
    for i = 2:length(u)-1
        d[i] = 0.5*(u[i+1] - u[i-1])/Δx
    end
    d[end] = 0.5*(u[1] - u[end-1])/Δx
    return d
end

"""
    transverse_step!(advance_fields!, em, pʸ, pᶻ, density, t, Δt)

One step of the transverse half of the reduced model: accumulate the canonical
`p⊥ = -A⊥` from `E⊥`, then advance `em` with the current that momentum carries.

Three lines, and every one of them has been wrong at some point, which is why
they are a function rather than a passage inside [`wakefield`](@ref):

* the accumulation is `p⊥ += E⊥Δt`, which is the invariant `p⊥ = -A⊥` and not a
  force integral -- see the `wakefield` docstring;
* the current argument is `-J Δt` rather than `J`, because
  `make_advance_fields` adds it straight into the field. Without the `Δt` the
  wakefield study reached a peak field of 1.0e22;
* the sign pairs `∂p/∂t = +e` with `∂e/∂t = -n·p` into an oscillation. The
  other way round it is exponential growth, measured at 44 rather than 1.

`density` is whatever the caller has: the Vlasov density of the step in
`wakefield`, a constant in `test_em_plasma.jl`. Sharing this function is what
lets the dispersion relation asserted there be a statement about the code the
wakefield study runs, rather than about a second copy of it.
"""
function transverse_step!(advance_fields!, em, pʸ, pᶻ, density, t, Δt)
    pʸ .= pʸ .+ em.ey.*Δt
    pᶻ .= pᶻ .+ em.ez.*Δt
    advance_fields!(t, (y = -pʸ.*density.*Δt, z = -pᶻ.*density.*Δt))
    return nothing
end

"""
    em_omega(k, density, Δx, Δt)

Frequency of a plane electromagnetic wave in uniform plasma **on this grid**:

    (2/Δt)²·sin²(ωΔt/2) = (2/Δx)²·sin²(kΔx/2) + n

Not an approximation and not the continuum `ω² = ωₚ² + k²` — it is exact for
the update order in [`transverse_step!`](@ref), which is worth setting out
because the two are far apart at the resolutions this package runs. Writing
`z = exp(-iωΔt)` and eliminating `H` and `p⊥` from the three update lines leaves

    -2i·sin(ωΔt/2)·E = -2i·cfl·sin(kΔx/2)·H - Δt²n·E/(-2i·sin(ωΔt/2))

and the relation above follows. The plasma enters only through `Δt²n/4`, so at
the small `cfl` the wakefield study uses it is the *spatial* term that carries
the error: at ten cells per wavelength the group velocity comes out 4.5% below
`√(1-n)`, which is a laser pulse arriving late and a wake with the wrong phase
velocity, not a rounding difference.

Real for `cfl²sin²(kΔx/2) + nΔt²/4 ≤ 1`, which is the stability condition; the
plasma term is what tightens it slightly below the familiar `cfl ≤ 1`.
"""
function em_omega(k, density, Δx, Δt)
    s = (Δt/Δx)^2*sin(k*Δx/2)^2 + density*Δt^2/4
    s ≤ 1 || error("em_omega: cfl²sin²(kΔx/2) + nΔt²/4 = $s exceeds 1, so this " *
                   "grid is unstable at k = $k and the frequency is complex")
    return 2/Δt*asin(sqrt(s))
end

"""
    vg_discrete(ω, density, Δx, Δt)

Group velocity `dω/dk` of [`em_omega`](@ref), differentiated rather than
estimated:

    v_g = (K/W)·cos(kΔx/2)/cos(ωΔt/2),    W = (2/Δt)sin(ωΔt/2), K = √(W² - n)

The continuum answer is `√(1 - n/ω²)`, and the two are not close on a coarse
grid: at ten cells per vacuum wavelength and `n = 0.1` this gives 0.906 against
0.949. Both factors matter and they pull the same way -- `cos(kΔx/2)` is the
Yee stencil running out of resolution, `K/W` is the plasma.
"""
function vg_discrete(ω, density, Δx, Δt)
    W = 2/Δt*sin(ω*Δt/2)
    W^2 > density || error("vg_discrete: ω = $ω is at or below the cutoff for n = $density")
    K = sqrt(W^2 - density)
    K*Δx/2 ≤ 1 || error("vg_discrete: ω = $ω does not propagate on a grid this coarse")
    return (K/W)*cos(asin(K*Δx/2))/cos(ω*Δt/2)
end

"""
    vg_pulse(density, Δx, Δt; duration, ω₀ = 1.0)

Group velocity of [`vg_discrete`](@ref) averaged over the spectrum of a pulse
with envelope `exp(-(ξ/duration)²)`, whose intensity spectrum is
`exp(-(ω - ω₀)²·duration²/2)`.

**The averaging is the small correction, and saying so is the point.** A pulse
one cycle long looks like it should have no single group velocity, and that
reading was once used to explain why the wakefield study's driver travelled at
0.886 against a continuum `√(1-n)` of 0.949. It does not: the same average over
the *continuum* `v_g` moves it from 0.9487 to 0.9438, half a percent. The other
six percent was the grid, which is what [`vg_discrete`](@ref) accounts for.

Measured against the study's own driver: 0.8993 predicted against 0.8858
measured at ten cells per wavelength, 0.9329 against 0.9261 at twenty.
"""
function vg_pulse(density, Δx, Δt; duration, ω₀ = 1.0)
    ωs = range(0.4*ω₀, 1.6*ω₀; length = 2001)
    w = @. exp(-(ωs - ω₀)^2*duration^2/2)
    return sum(w .* vg_discrete.(ωs, density, Δx, Δt))/sum(w)
end

"""
    wakefield(; Δx, Δt_factor, Δp, total_time, ...)

Laser wakefield excitation in a 1D1V plasma slab: `PFC` advection in x and p,
a spectral Poisson solve for the longitudinal field, `FDTD1D` for the
transverse one, and the ponderomotive force that couples the two.

Returns `(; t, x, p, nᵢ, n, ey, ex, Φ, ε_e, ε, plasma_density,
plasma_temperature)` -- the ion background, the density, transverse and
longitudinal field and ponderomotive-potential histories as `Nt × Nx` matrices,
the electrostatic and total energy histories, and the two plasma parameters the
theory the test compares against is built from.

# How the laser drives the wake

The transverse momentum is the **canonical** one rather than a force integral.
In this geometry nothing depends on `y` or `z`, so `p⊥ + A⊥` is conserved along
a trajectory; the plasma is at rest ahead of the pulse, where `A⊥` is zero, so
the constant is zero and `p⊥ = -A⊥` everywhere the pulse has reached. With
`E⊥ = -∂A⊥/∂t` that makes the accumulation below, `p⊥ += E⊥Δt`, the local
vector potential and not an equation of motion -- the same line either way, and
worth saying, because read as a force integral it would be missing the `v×B`
term and the convective derivative, and read as the invariant it is exact.

The longitudinal push then carries the ponderomotive force,

    ∂p/∂t = e - ∂Φ/∂x,       Φ = (pʸ² + pᶻ²)/2

and that term is what makes this a wakefield run at all. Without it the laser
never entered the longitudinal dynamics and `ex` was the slab edges relaxing --
`peak wake field` and `Δε/ε` came out bit-identical whether the transverse
current was right, wrong by `Δt`, or wrong by thirty-two orders of magnitude,
which is what this function computed until the term was added.

`Φ` is the non-relativistic potential: the relativistic `(pʸ² + pᶻ²)/2γ` varies
along `p`, and the momentum sweep displaces a whole column by one scalar, so
carrying `1/γ` would mean a different advection rather than a different number.
Both that and the current, which is taken through momentum rather than velocity,
are corrections of relative order `p²`, and the verification runs at an
amplitude where that is a fraction of a percent. Neither is a licence to read
this as a relativistic solver.

# The defaults are the resonance, and the amplitude is bounded by the model

The wake follows the spectrum of `Φ` at the wake number `k_p`, so a pulse much
longer than `1/k_p` drives nothing: the Gaussian factor is `exp(-k_p²σ²/2)`, and
the `laser_duration = 5·2π` this study ran at before, at `k_pσ = 5.7`,
suppressed the wake by seven orders of magnitude. That is the second reason the
wake it produced was entirely the slab edges -- there was no term to drive one,
and the pulse was the wrong length to drive one had there been.
`laser_duration = 2π` sits at `k_pσ = 1.13`, measuring `k_p` as the wake comes
out at `0.360`.

`laser_amplitude = 0.3` is likewise a statement about the model rather than a
taste: it reaches `p⊥ = 0.26`, where the `1/γ` this drops is worth 3.2%, and a
wake 2.8% of the cold wave-breaking field `√n`. At the `1.0` that stood here
before those are 25% and 31%, and neither the ponderomotive potential nor the
current is the right expression any more.

# `Δx` resolves the *laser*, and twenty cells is the floor

`Δx = 0.05·2π` is twenty cells per vacuum wavelength. It was ten, and ten is not
a resolution this model can be read at: the Yee dispersion relation
[`em_omega`](@ref) puts the driver's group velocity 4.5% below `√(1 - n)` there,
so the pulse arrives late and the wake it writes has a phase velocity wrong by
the same amount -- `γ_φ ≈ 2.2` against the 3.0 the physics gives, which is the
difference between two different statements about trapping and dephasing.
Measured across the two, at about four times the wall clock -- halving `Δx` halves
`Δt` with it, so `Nx` and `Nt` both double (2.5 s against 10 s, timed back to back):

    cells/λ₀   pulse speed   wake phase velocity   λ       predicted v_g
    10         0.8858        0.8925                17.460  0.8993
    20         0.9261        0.9207                18.010  0.9329

The wake's *frequency* is insensitive to this -- 19.564 against 19.561 for the
period -- which is what makes the error easy to miss: every assertion about the
oscillation passes at either resolution, and only the ones about the wavelength
and the pulse speed move. Nothing here is converged in `Δx` to better than a few
percent even now; `test_verification.jl` asserts the trend rather than pretending
otherwise.

Otherwise the defaults are the study's own parameters -- the script plots what
they produce, the test asserts it -- and they are not to be coarsened for speed:
a coarsened run is a different experiment, not a faster version of this one.
"""
function wakefield(; Δx = 0.05*2π,
                     Δt_factor = 0.05,
                     Δp = 0.02,
                     total_time = 2π*22,
                     x_min = -5.0*2π,
                     box_length = 20.0*2π,
                     plasma_thickness = 10.0*2π,
                     plasma_temperature = 0.01,
                     plasma_density = 0.1,
                     laser_amplitude = 0.3,
                     laser_duration = 2π,
                     p_max = 5*sqrt(plasma_temperature))
    Δt = Δt_factor*Δx
    x = collect(x_min:Δx:(box_length + x_min))
    p = collect(-p_max:Δp:p_max)
    Nx, Np = length(x), length(p)

    f = plasma_density/sqrt(2π*plasma_temperature) *
        (@. exp(-0.5*(p)^2/plasma_temperature)) *
        (@. 0.5*(tanh(x) - tanh(x - plasma_thickness)))'

    nᵢ = integrate(p, f)          # immobile neutralising ions
    g = similar(f)

    # PFC holds no arrays, so one value serves every line of both sweeps.
    advection = PFC(fmin = 0.0, fmax = maximum(f))

    em = FDTD1D.YeeMesh1D{Float64}(Nx - 1)
    # The pulse is launched by its own peak rather than by its tail. The source
    # is driven at `x_min + Δx` with a profile in `x - t`, so writing it in
    # `x - t` alone puts the maximum at `t = x_min + Δx`, which is *before* the
    # run starts: this study injected `exp(-((x-t)/L)²)` from `t = 0` at
    # `x_min = -5·2π`, reached 0.383 of its own amplitude, and called that the
    # peak laser field. The delay moves the maximum to `t = laser_delay`, and
    # three durations of it leave `exp(-9)` on the grid at `t = 0`.
    laser_delay = 3*laser_duration
    ξ(t, x) = x - t - (x_min + Δx - laser_delay)
    pulse_shape = (y = (t, x) -> laser_amplitude*exp(-(ξ(t, x)/laser_duration)^2)*sin(ξ(t, x)),
                   z = (t, x) -> 0.0)
    advance_fields! = FDTD1D.make_advance_fields(
        em, Δt/Δx, pulse_shape, Δt, Δx, x_min,
        FDTD1D.PML(; N = 0, σ_max = 1.0, Δx = Δx, Δt = Δt))

    t = collect(0.0:Δt:total_time)
    Nt = length(t)
    n = zeros(Nt, Nx)
    ey = zeros(Nt, Nx)
    ex = zeros(Nt, Nx)
    Φ = zeros(Nt, Nx)
    ε_e = zeros(Nt)
    ε = zeros(Nt)

    solve_poisson! = PoissonFourier1D.generate_solver(nᵢ, Δx)
    e = similar(x)
    n[1, :] = nᵢ
    solve_poisson!(e, n[1, :] - nᵢ)
    ex[1, :] = e
    ey[1, :] = em.ey
    ε_e[1] = integrate(x, e.^2)
    ε[1] = integrate(x, integrate(p, @. f*p^2)) + ε_e[1]

    pʸ = zeros(Nx)
    pᶻ = zeros(Nx)
    ϕ = zeros(Nx)                 # the ponderomotive potential of this step
    ∂ϕ = zeros(Nx)
    force = zeros(Nx)

    for k in 2:Nt
        for j = 1:Np
            advect!(view(g, j, :), view(f, j, :), advection, p[j]*Δt/Δx)
        end

        n[k, :] = integrate(p, g)
        solve_poisson!(e, n[k, :] - nᵢ)

        # The canonical transverse momentum and the current it carries, both in
        # `transverse_step!` so that `test_em_plasma.jl` can assert the
        # dispersion relation of *this* code rather than of a copy of it.
        transverse_step!(advance_fields!, em, pʸ, pᶻ, view(n, k, :), k*Δt, Δt)

        @. ϕ = 0.5*(pʸ^2 + pᶻ^2)
        ddx!(∂ϕ, ϕ, Δx)
        # The ponderomotive force expels from high intensity whatever the sign
        # of the charge, so it enters against the gradient while `e` enters
        # with it -- the one term here whose sign is fixed by the physics
        # rather than by this package's field convention.
        @. force = e - ∂ϕ

        for i = 1:Nx
            advect!(view(f, :, i), view(g, :, i), advection, force[i]*Δt/Δp)
        end

        ey[k, :] = em.ey
        ex[k, :] = e
        Φ[k, :] = ϕ
        ε_e[k] = integrate(x, e.^2)
        ε[k] = integrate(x, integrate(p, @. f*p^2)) + ε_e[k]
    end

    # The two plasma parameters come back with the run. The test needs them to
    # form `ωₚ` and the Bohm--Gross shift, and reading them off the result is
    # what keeps it measuring the run it was handed rather than the defaults it
    # was written against.
    return (; t, x, p, nᵢ, n, ey, ex, Φ, ε_e, ε, plasma_density, plasma_temperature)
end

"""
    linear_wake(x, t, Φ, nᵢ; temperature)

The wake linear theory puts behind the recorded ponderomotive drive `Φ`.

This is the closed form the test measures `wakefield` against, and it shares no
code with it. `Φ` is the laser's, off the `FDTD1D` side; the wake it is compared
with comes off the Vlasov push and the Poisson solve; this ODE is the line
between them, and it is the only thing claiming they should agree. Linearising
the warm fluid about a fixed ion background, in this package's conventions --
`∂p/∂t = e` for the push, `∂e/∂x = n - nᵢ` for what `PoissonFourier1D` returns,
and `P₁ = 3Tn₁` for the one-dimensional adiabatic closure --

    ∂n₁/∂t + n₀ ∂u/∂x = 0,   ∂u/∂t = e - ∂Φ/∂x - 3T ∂n₁/∂x / n₀,   ∂e/∂x = n₁

and eliminating `n₁` and `u` leaves a plasma oscillator driven by the laser:

    ∂²e/∂t² + ωₚ² e - 3T ∂²e/∂x² = ωₚ² ∂Φ/∂x,      ωₚ²(x) = nᵢ(x)

**Driven by the recorded `Φ(x,t)` rather than by an idealised pulse**, which is
what lets the comparison be pointwise. A reference built from the pulse's
analytic envelope would have to assume the shape the laser arrives with, that it
translates rigidly, and that it neither disperses nor gives up a measurable
fraction of itself to the wake; taking the drive from the run under test assumes
none of that, and what remains to be verified -- whether the plasma responds to
that drive as theory says -- is the part `wakefield` computes.

**The `3T` term is not decoration.** It is the Bohm--Gross correction that
`test/test_verification.jl` already measures against `√(1+3k²)` for a standing
oscillation, and dropping it here leaves the reference oscillating 1.7% slow.
That is a third of a radian over the three plasma periods a point in the window
has been ringing for by the end of the run, and it takes the agreement from 8.9% to
26% -- the comparison fails on phase long before anything is wrong with the
amplitude, which the same change moves by 6%. Near the slab edge the term
carries a `∂n₀/∂x` this drops, so the reference is a statement about the
interior.
"""
function linear_wake(x, t, Φ, nᵢ; temperature)
    Δx = x[2] - x[1]
    e = zeros(length(x))
    u = zeros(length(x))          # ∂e/∂t
    ∂Φ = similar(e)
    ∂e = similar(e)               # `ddx!` reads its source per point, so the
    ∂²e = similar(e)              # two passes cannot share one buffer
    for k in 2:length(t)
        Δt = t[k] - t[k-1]
        ddx!(∂Φ, view(Φ, k, :), Δx)
        ddx!(∂e, e, Δx)
        ddx!(∂²e, ∂e, Δx)
        @. u += Δt*(nᵢ*(∂Φ - e) + 3*temperature*∂²e)
        @. e += Δt*u
    end
    return e
end

"""
    wake_wavelength(v_drive, density, temperature)

The wake's spatial period behind a driver moving at `v_drive`.

The wave that stays in step with the driver is the one whose phase velocity
*is* the driver's, so `ω(k) = k·v_drive` against the Bohm--Gross
`ω² = ωₚ² + 3Tk²` fixes `k`, and

    λ = 2π√(v_drive² - 3T)/ωₚ

The thermal term shortens it by 1.7% at the defaults here, and enters with the
opposite sign to the way it enters a standing oscillation: there it raises the
frequency at fixed `k`, here the frequency is pinned by the driver and it raises
`k` instead.

**`v_drive` is the driver's group velocity on the grid, not `√(1 - n)`.** This
docstring used to say there was no closed form to reach for: the continuum
`√(1 - n)` is 0.949 at this density while the pulse travelled at 0.886, and the
gap was put down to a one-cycle driver having a bandwidth of order its own
carrier. That explanation was wrong, and measurably so. Averaging the continuum
group velocity over the pulse's own spectrum moves it from 0.9487 to 0.9438 --
half a percent, not six. The rest was the Yee grid: at the ten cells per
wavelength that study ran, [`vg_pulse`](@ref) gives 0.8993 against the 0.8858
measured, and halving `Δx` moves the measurement to 0.9261 against a predicted
0.9329. A missing six percent in the driver is a wake whose phase velocity is
wrong by six percent, which is `γ_φ ≈ 2.2` where the physics says 3.0.

So the honest division is not "measure the driver, predict the wake" but
"predict both": [`vg_pulse`](@ref) is a closed form for the driver, this is one
for the wake it writes, and [`pulse_velocity`](@ref) measures the first
independently so the two can be compared rather than assumed.
"""
wake_wavelength(v_drive, density, temperature) =
    2π*sqrt(v_drive^2 - 3*temperature)/sqrt(density)

"""
    pulse_velocity(t, x, Φ; lo, hi)

Speed of the ponderomotive peak, by least squares over the samples where it lies
in `lo ≤ x ≤ hi`.

The window is there to fit the crossing of the slab interior only: the peak
tracked outside it is a pulse still entering the plasma, or one already leaving,
and neither travels at the speed the wake behind it was written at.
"""
function pulse_velocity(t, x, Φ; lo, hi)
    xs, ts = Float64[], Float64[]
    cut = 0.2*maximum(Φ)
    for k in 2:length(t)
        i = argmax(view(Φ, k, :))
        if lo ≤ x[i] ≤ hi && Φ[k, i] > cut
            push!(xs, x[i]); push!(ts, t[k])
        end
    end
    length(ts) ≥ 10 || error("pulse_velocity: the peak spent $(length(ts)) samples in [$lo, $hi]")
    return (hcat(ones(length(ts)), ts) \ xs)[2], length(ts)
end

"""
    wave_period(u, y; lo, hi)

Period of `y` against `u`, from the mean spacing of its interpolated zeros,
doubled -- zeros come twice a period.

`u` is whichever axis is being measured along: pass `x` and it returns the
wake's wavelength, pass `t` and it returns the wake's oscillation period. The
two are the same measurement, and keeping them one function is what makes the
phase velocity the ratio of two things measured the same way rather than two
conventions that have to be reconciled.

Zeros rather than the extrema [`local_extrema`](@ref) finds, because an extremum
is a grid point and a spacing built from grid points is quantised to `Δx` --
3.3% of this wake's period at the study's resolution, which is coarser than the
1.7% thermal shift being resolved.
"""
function wave_period(u, y; lo, hi)
    z = zero_crossings(u, y; lo = lo, hi = hi)
    length(z) ≥ 3 || error("wave_period: $(length(z)) zeros in [$lo, $hi], need 3")
    return 2*(z[end] - z[1])/(length(z) - 1), length(z)
end

"""
    zero_crossings(x, y; lo, hi)

Interpolated zeros of `y` in `lo ≤ x ≤ hi`, for measuring a wavelength to better
than the grid.

The extrema `local_extrema` returns are grid points, and a spacing built from
them is quantised to `Δx` -- 3.3% of the wake's period at this study's
resolution, which is coarser than the 1.7% thermal shift the wavelength is
supposed to resolve. A zero is a crossing between two samples and interpolates
linearly to a fraction of a cell.
"""
function zero_crossings(x, y; lo, hi)
    z = Float64[]
    for i in 1:length(y)-1
        (lo ≤ x[i] && x[i+1] ≤ hi) || continue
        (y[i] == 0 || sign(y[i]) != sign(y[i+1])) || continue
        push!(z, x[i] - y[i]*(x[i+1] - x[i])/(y[i+1] - y[i]))
    end
    return z
end

"""
    growth_rate(t, ε_e; lo, hi)

Instability growth rate `γ` from `ε_e ∝ exp(2γt)`, fitted by least squares over
the stretch where `ε_e` rises from `lo` to `hi`.

**A plain fit over every sample, unlike [`damping_rate`](@ref).** That is not an
inconsistency: the unstable root of the cold two-stream dispersion relation is
*purely imaginary*, so the mode grows without oscillating and `log ε_e` is a
straight line with no `log cos²` poles to fall into. The Landau mode has a real
frequency and needs its envelope; this one has none and does not. There are in
fact no local maxima here to fit through, so `damping_rate` would raise on the
first call.

**The window is set by amplitude rather than by time**, which is what makes it
transferable between wavenumbers: the growth rate varies over the branch, so a
fixed time window covers a different stretch of the exponential at each `k` and
the fitted value wobbles by several percent with it. Measured at `kv₀ = 0.4`
over the same run, fitting `t ∈ [8,18]`, `[10,20]`, `[12,22]`, `[14,24]` gives
9.75%, 5.63%, 4.35% and 0.55% error; the amplitude band gives 1.86% and does the
same thing at every `k`.

!!! note "Why a fixed window wobbles: `ε_e` is not one exponential"
    The quadratic in `γ_cold` has **four** roots -- the growing pair `±iγ` and
    an oscillating pair `±ω₊`, with `ω₊ = √u₊` and
    `u₊ = [(2a² + 1) + √(8a² + 1)]/2`. An initial perturbation excites all of
    them, and the cross term between the growing root and the oscillating ones
    puts a ripple on `ε_e` at `ω₊` whose size **relative to** the growing mode
    falls only as `exp(−γt)`. So it is still there through any window one can
    afford to fit over.

    Measured at `a = 0.6`: the instantaneous rate oscillates with period 4.5
    against the `2π/ω₊ = 4.626` this predicts, swinging between 0.21 and 0.41
    around a `γ_cold` of 0.353. Fitting over an integer number of beat periods
    instead of an arbitrary window cuts the spread over start points from
    39.6%, 14.4% and 22.3% (at `a` = 0.4, 0.6, 0.8) to 9.0%, 5.6% and 4.8%.

    The ripple is worst where `γ` is smallest, since that is what sets how fast
    it decays away -- which is why `a = 0.4` and `a = 0.9`, at either end of the
    branch, scatter more than `a = 0.6` and `0.8` near the peak.

    **The amplitude band already handles this**, which is the reason not to do
    anything cleverer: it spans 1.15 to 2.22 beat periods across the three cases
    in use, enough to average the ripple. Adding `cos ω₊t` and `sin ω₊t` to the
    design matrix -- still a linear fit, since `ω₊` is known in closed form --
    was tried and moves the band results by at most one point (−1.86% to
    −2.86%, −3.14% to −3.08%, +0.31% to +0.44%). It is not worth the machinery.

    This is what produced the apparent overshoot above `γ_cold` at small beam
    temperature: `vt` changes `γ` slightly, which moves the beat's phase within
    a fixed window, and the fitted rate follows it across the cold value. Two
    other explanations were measured and rejected first -- refining `Δv` moves
    the result by 1e-5, and the driver's renormalisation leaves the effective
    density at 1.0000158, worth 0.0008% on `γ`. So was a third: at `a = 0.4` the
    second harmonic really is more unstable than the fundamental
    (`γ(0.8) = 0.311` against `γ(0.4) = 0.308`), but it starts at `O(α²)` and
    gains 13% over the run against a head start of 1e-6, so it contributes
    nothing here.

`hi` also has to keep the run inside the solver's validity. The field grows with
the mode, and the velocity sweep is displaced by `E·Δt`, so a large enough `ε_e`
breaks `PFC`'s Courant limit in `v` and the run diverges -- measured `ε_e` at
1.2e161 before `NaN` at `t = 24.1`. Widening the velocity window only postpones
it, from `t = 24.1` at `±8` to `t = 26.6` at `±16`, which is what identifies the
Courant limit rather than the boundary as the cause. At the `hi = 5.0` this
package uses, the velocity Courant number `max|e|Δt/Δv` is 0.52 to 0.66 in the
three growth-rate runs and 0.73 at `a = 1.0`. (It read 0.46 to 0.65: that is a
single mode's amplitude `√(2ε_e/L)`, and the peak of `e` is higher.)
[`two_stream`](@ref) stops its runs at 1.
"""
function growth_rate(t, ε_e; lo, hi)
    i0 = findfirst(≥(lo), ε_e)
    i1 = findfirst(≥(hi), ε_e)
    (i0 === nothing || i1 === nothing) &&
        error("growth_rate: ε_e never spans [$lo, $hi] (range $(extrema(ε_e)))")
    i1 - i0 ≥ 10 ||
        error("growth_rate: only $(i1 - i0 + 1) samples between $lo and $hi")
    band = @view ε_e[i0:i1]
    # A run that breaches the Courant limit early enough puts a non-finite
    # sample *inside* the band rather than after it, and `A \ log.(...)` then
    # returns a `NaN` slope that fails a downstream `isapprox` with nothing to
    # point at. Diagnose it here, where the cause is still visible. `≤ 0` is
    # caught with it: `log` of a zero sample would give `-Inf` and the same
    # silent `NaN`.
    all(x -> isfinite(x) && x > 0, band) || error(
        "growth_rate: the fit window t ∈ [$(t[i0]), $(t[i1])] contains a " *
        "non-positive or non-finite ε_e (first at t = " *
        "$(t[i0 + findfirst(x -> !(isfinite(x) && x > 0), band) - 1])). " *
        "The run has diverged into the band being fitted -- shorten it, or " *
        "lower `hi` so the fit ends before the velocity Courant limit.")
    A = hcat(ones(i1 - i0 + 1), t[i0:i1])
    return (A \ log.(band))[2]/2, t[i0], t[i1]
end

# -------------------------------------------------- two-stream instability
#
# Shared rather than duplicated, for the reason `wakefield` is: the plotting
# script and the test that asserts its claims have to run the same setup. The
# `tmax` note below is the whole argument -- a copy of this function with the
# number and without the paragraph explaining why it cannot move is a trap, and
# that is what `verification/two-stream.jl` carried until this was hoisted.

"""
    γ_cold(a)

Growth rate of the **cold** two-stream instability at `a = kv₀`, for two beams
of density 1/2 at `±v₀` with `ω_p = 1`.

Derived rather than quoted. The electrostatic dispersion relation is

    1 = ½/(ω − kv₀)² + ½/(ω + kv₀)²

which, with `u = ω²`, is a quadratic:

    u² − (2a² + 1)u + (a⁴ − a²) = 0
    u± = [(2a² + 1) ± √(8a² + 1)]/2

`u₋ < 0` exactly when `a < 1`, and then `γ = √(−u₋)`. So this case needs no
plasma dispersion function, no numerical root and no tabulated constant of the
kind the Landau cases have to carry -- which is what made it worth waiting for
rather than hard-coding a number. `test_verification.jl` checks the closed form
against the relation it came from before using it.

**This is the `vt → 0` limit, not the case the runs are held to.** The beams in
[`two_stream`](@ref) are Maxwellian at `vt = 0.3`, where the cold rate is off by
up to 3.14% -- and by 7.44% for the `vt = 0.6` run, the cold error growing with
temperature; [`two_stream_warm`](@ref) solves the same relation for warm beams
and is what the assertions compare against. The two meet to 0.006% at
`vt = 0.02`, which `test_dispersion.jl` asserts, and they part company in sign
as well as size near the band edge: above `a ≈ 0.77` the warm rate is the larger
of the two, and above `a = 1` the cold branch is zero while warm beams still
grow.
"""
two_stream_u(a) = ((2a^2 + 1) - sqrt(8a^2 + 1))/2
γ_cold(a) = sqrt(max(0.0, -two_stream_u(a)))

"The cold dispersion relation itself, for checking `γ_cold` against."
two_stream_residual(ω, a) = 0.5/(ω - a)^2 + 0.5/(ω + a)^2 - 1

"""
    two_stream(a; v₀, vt, Δv, vmax, Δt, tmax)

Two counter-streaming warm beams at `±v₀`, perturbed by 0.1% in the `k = a/v₀`
mode, returning `(t, ε_e)` over one wavelength.

Top level rather than a closure inside the testset, matching `landau_case`.
That is a readability choice and not a performance one: moving it out was tried
as a fix for what looked like 50 s of compilation, and changed the runtime not
at all. The 50 s was a mismeasurement -- the block costs 6.7 s, of which 3.8 s
is arithmetic, against a testset that already took 1m12 before it was added.

`vmax = 6` rather than the 8 first used. The growth rate is a linear-phase
measurement, taken before the beams have spread, so it does not see the window
at all -- measured identical to five digits at `vmax` 5, 6 and 8 -- and the
narrower grid halves the cost.

!!! note "The runs stop at the velocity Courant limit, and `tmax = 24.0` has only to be long enough"
    An unstable run grows its field until the velocity sweep moves `f` by more
    than a cell per step, and past that `PFC` leaves the bounds it is built on.
    Unchecked, the `a = 0.6` run crossed the limit at `t = 21.1`, handed the
    scheme `f = −5.5e-20` at `t = 22.0`, with the Courant number at 1.41, and
    −3.9e-3 by 23.4, and went non-finite at 24.15. The scheme refuses that now,
    so `two_stream` runs with `stop_at_courant = true` and each run ends at its
    own limit: `a = 0.4` with the step from `t = 23.2`, `0.6` from 21.1, and the
    `vt = 0.6` run from 22.2. `a = 0.8` does not reach it by `t = 24` (0.94), nor
    `a = 1.0` by 80 (0.94).

    The fits are over long before: `hi = 5` holds the Courant number at 0.52 to
    0.66 in these runs, and at 0.73 in the `a = 1.0` one, and none of the rates
    moved when the runs started stopping. What bounds `tmax` is the slowest fit
    completing -- the `a = 0.8` one needs `ε_e` to reach 5.0, which happens at
    `t = 22.85`, and below that `growth_rate` raises rather than guessing. It
    used to be bounded above as well, by the fastest run diverging, which left a
    window of three steps; that bound is gone. `growth_rate` still raises on a
    window containing a non-finite sample, for a run set up to cross the limit
    inside its own fit.

    **`a = 1.0` is a different regime and takes `tmax = 80`.** There the cold
    rate is zero and the warm one is 0.098, a third of the branch maximum, so
    `ε_e` needs 78 time units to cross the same amplitude band the other cases
    cross in twenty.
"""
function two_stream(a; v₀ = 3.0, vt = 0.3, Δv = 0.05, vmax = 6.0,
                       Δt = 0.05, tmax = 24.0)
    k = a/v₀
    L = 2π/k
    Nx = round(Int, L/0.49)
    Δx = L/Nx
    x = collect(Δx:Δx:L)
    v = collect(-vmax:Δv:vmax)
    t = collect(0.0:Δt:tmax)
    beams = @. 0.5/sqrt(2π*vt^2)*(exp(-(v - v₀)^2/(2vt^2)) +
                                  exp(-(v + v₀)^2/(2vt^2)))
    f₀ = beams * (@. (1.0 + 1e-3*cos(k*x)))'
    r = vlasov_poisson(x, v, f₀, t; stop_at_courant = true)
    n = length(r.ε_e) - 1
    return t[1:n], r.ε_e[1:n]
end

