# Shared 1D1V Vlasov–Poisson driver for the extended verification tests.
# Deliberately close to what the verification notebooks do, so the two agree.

using Vasilek
using Vasilek: StrangSplitting, FDTD1D, PoissonFourier1D
using NumericalIntegration, FFTW

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
    vlasov_poisson(x, v, f₀, t; scheme_x, scheme_v, invariants = false)

Strang-split Vlasov–Poisson on a static grid.

Returns a NamedTuple with `ε_e` (electric energy) and `ε` (total energy)
histories. With `invariants = true` it also returns the `mass`, `momentum`,
`l2` and `entropy` histories -- the conserved quantities that are *not* the
energy, and that nothing asserted until now.

**The invariants use the cell-width sum `Σ f ΔvΔx`, not `integrate`.** That is
the quadrature the schemes actually conserve: `PFC` is a flux form, so what
leaves one cell enters its neighbour and the full-weight sum is preserved
exactly. The trapezoid halves the two endpoint weights, which no flux
conservation law protects, and measuring with it reports a drift that belongs to
the quadrature rather than to the scheme. Measured over the k = 0.5 Landau run,
875 steps: mass drifts 2.8e-16 by the cell-width sum against **1.6e-4** by the
trapezoid, and momentum stays at 1.7e-15 against 7.3e-4. Both trapezoid figures
are the endpoint weighting, not the solver.

The energy histories above keep `integrate`, because they are compared with
tolerances of half a percent where the difference is irrelevant, and because
changing them would silently move numbers the notebooks quote.

`scheme_x` and `scheme_v` default to `PFCNonUniform` on the two grids, which is
what the verification notebooks use and what every previous caller got. They are
arguments so that the same driver can measure what the physics costs under a
*different* scheme, which is what `verification/scheme-comparison.jl` does, and
so that a refinement study can hold the scheme fixed while moving the grid.
"""
function vlasov_poisson(x, v, f₀, t;
                        scheme_x = nothing, scheme_v = nothing, invariants = false)
    Δx = cell_widths(x)
    Δv = cell_widths(v)
    sx = scheme_x === nothing ? PFCNonUniform(Δx; fmin = 0.0, fmax = 1.0) : scheme_x
    sv = scheme_v === nothing ? PFCNonUniform(Δv; fmin = 0.0, fmax = 1.0) : scheme_v

    advect_x! = line_advector(sx, Δx)
    advect_v! = line_advector(sv, Δv)

    fᵢ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. Δx/Δx)'
    nᵢ = integrate(v, fᵢ)
    Nᵢ = integrate(x, nᵢ)

    f = copy(f₀)
    f .*= Nᵢ/integrate(x, integrate(v, f))
    g = f'

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

    for k in 1:length(t)-1
        Δt = t[k+1] - t[k]
        vΔt(_) = v*Δt
        function eΔt(ff)
            nₖ = vec(sum(ff'.*Δv, dims = 1))
            solve_poisson!(e, nₖ - nᵢ)
            return e*Δt
        end
        StrangSplitting.make_time_step_2d!((g, f), (vΔt, eΔt), (advect_x!, advect_v!))
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
        end
    end
    for h in (ε_e, ε, mass, momentum, l2, entropy, fmin)
        h === nothing || (h[end] = h[end-1])
    end
    return (; ε_e, ε, mass, momentum, l2, entropy, fmin)
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

Otherwise the defaults are the study's own parameters -- the script plots what
they produce, the test asserts it -- and they are not to be coarsened for speed:
a coarsened run is a different experiment, not a faster version of this one.
"""
function wakefield(; Δx = 0.1*2π,
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

        # The canonical transverse momentum, `p⊥ = -A⊥`. See the docstring.
        pʸ .= pʸ .+ em.ey.*Δt
        pᶻ .= pᶻ .+ em.ez.*Δt

        # The current owes its own Δt -- `make_advance_fields` adds the argument
        # straight into `ey` -- and the sign that pairs `∂p/∂t = +e` with
        # `∂e/∂t = −n·p` into an oscillation rather than exponential growth.
        # Without the first the peak field reached 1.0e22; with the sign the
        # other way, 44. See `docs/normalization.md`.
        advance_fields!(k*Δt, (y = -pʸ.*n[k, :].*Δt, z = -pᶻ.*n[k, :].*Δt))

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

**`v_drive` is measured rather than assumed, and there is no closed form to
reach for instead.** The obvious candidate is the monochromatic group velocity
`v_g = √(1 - n)`, which is 0.949 at this density; the pulse actually travels at
0.886, because a driver one cycle long has a bandwidth of order its own carrier
and no single `ω₀` describes it. That gap is not absorbable at these
tolerances, and the arithmetic is worth writing down because `√(1-n)` is a
tempting thing to substitute: it gives `λ = 18.533` against a measured 17.460,
an error of 6.14% against the `rtol = 0.03` the test holds `λ` to -- more than
double the tolerance, so the assertion fails rather than drifting. A
`group_velocity` helper computing `√(1-n)` used to sit above this function,
unused and recommending exactly that substitution; it was deleted rather than
documented, there being no caller it could serve.

Taking `v_drive` from [`pulse_velocity`](@ref) leaves the closed form above as
the claim and the driver as an input to it -- which is the honest division, the
wake being the part this package computes.
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
9.76%, 5.63%, 4.35% and 0.56% error; the amplitude band gives 1.87% and does the
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
    was tried and moves the band results by at most one point (−1.87% to
    −2.87%, −3.14% to −3.08%, +0.31% to +0.44%). It is not worth the machinery.

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
package uses, the velocity Courant number is 0.46 to 0.65.
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
