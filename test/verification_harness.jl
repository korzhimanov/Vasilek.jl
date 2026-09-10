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
    wakefield(; Δx, Δt_factor, Δp, total_time, ...)

Laser wakefield excitation in a 1D1V plasma slab: `PFC` advection in x and p,
a spectral Poisson solve for the longitudinal field, and `FDTD1D` for the
transverse one.

Returns `(; t, x, p, n, ey, ex, ε_e, ε)` -- the density, transverse and
longitudinal field histories as `Nt × Nx` matrices, and the electrostatic and
total energy histories.

The defaults reproduce the study exactly, and every caller here uses them: the
script plots what they produce and the test asserts it. The keywords exist to
name the study's parameters and make them reachable -- for a resolution sweep,
or for the ponderomotive work the warning below describes -- and **not** so that
the test can run something cheaper.

That distinction is the reason the test takes a second rather than a tenth of
one. A coarsened run is a different experiment, not a faster version of this
one: doubling `Δx` alone moves the peak laser field from 0.383 to 0.641, and
doubling `Δx` and `Δp` together gives 0.641 with a 6.95% energy drift against
1.19%. The peak laser field is the single number that sees the transverse
current at all, so a proxy that moves it by two thirds would pin its own value
and call it the study's -- which is precisely the assertion the test exists to
make. Coarsen for exploration; do not coarsen and then assert.

!!! warning "Known incomplete"
    There is **no ponderomotive coupling**. The laser never enters the
    longitudinal push, so the wake this produces is the slab edges relaxing
    rather than a laser-driven wave -- `peak wake field` and `Δε/ε` come out
    bit-identical whether the transverse current is right, wrong by `Δt`, or
    wrong by thirty-two orders of magnitude. Closing that needs the
    ponderomotive force `−∇(pʸ² + pᶻ²)/2γ` in the momentum advection, which is a
    modelling decision rather than a repair.

    Anything asserted about the output is therefore a statement that the solver
    runs and stays bounded, not that the physics is complete. The test says so
    too, so that a passing run is not read as more than it is.

    Note also that the current is taken through momentum rather than velocity.
    At `laser_amplitude = 1.0` the motion is relativistic and `1/γ` is not close
    to unity, with `γ = sqrt(1 + pₓ² + pʸ² + pᶻ²)`.
"""
function wakefield(; Δx = 0.1*2π,
                     Δt_factor = 0.05,
                     Δp = 0.1,
                     total_time = 2π*10,
                     x_min = -5.0*2π,
                     box_length = 20.0*2π,
                     plasma_thickness = 10.0*2π,
                     plasma_temperature = 0.1,
                     plasma_density = 0.1,
                     laser_amplitude = 1.0,
                     laser_duration = 5*2π)
    Δt = Δt_factor*Δx
    x = collect(x_min:Δx:(box_length + x_min))
    p = collect(-laser_amplitude*4:Δp:laser_amplitude*4)
    Nx, Np = length(x), length(p)

    f = plasma_density/sqrt(2π*plasma_temperature) *
        (@. exp(-0.5*(p)^2/plasma_temperature)) *
        (@. 0.5*(tanh(x) - tanh(x - plasma_thickness)))'

    nᵢ = integrate(p, f)          # immobile neutralising ions
    g = similar(f)

    # PFC holds no arrays, so one value serves every line of both sweeps.
    advection = PFC(fmin = 0.0, fmax = maximum(f))

    em = FDTD1D.YeeMesh1D{Float64}(Nx - 1)
    pulse_shape = (y = (t, x) -> exp(-((x - t)/laser_duration)^2)*sin(x - t),
                   z = (t, x) -> 0.0)
    advance_fields! = FDTD1D.make_advance_fields(
        em, Δt/Δx, pulse_shape, Δt, Δx, x_min,
        FDTD1D.PML(; N = 0, σ_max = 1.0, Δx = Δx, Δt = Δt))

    t = collect(0.0:Δt:total_time)
    Nt = length(t)
    n = zeros(Nt, Nx)
    ey = zeros(Nt, Nx)
    ex = zeros(Nt, Nx)
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

    for k in 2:Nt
        for j = 1:Np
            advect!(view(g, j, :), view(f, j, :), advection, p[j]*Δt/Δx)
        end

        n[k, :] = integrate(p, g)
        solve_poisson!(e, n[k, :] - nᵢ)

        pʸ .= pʸ .+ em.ey.*Δt
        pᶻ .= pᶻ .+ em.ez.*Δt

        # The current owes its own Δt -- `make_advance_fields` adds the argument
        # straight into `ey` -- and the sign that pairs `∂p/∂t = +e` with
        # `∂e/∂t = −n·p` into an oscillation rather than exponential growth.
        # Without the first the peak field reached 1.0e22; with the sign the
        # other way, 44. See `docs/normalization.md`.
        advance_fields!(k*Δt, (y = -pʸ.*n[k, :].*Δt, z = -pᶻ.*n[k, :].*Δt))

        for i = 1:Nx
            advect!(view(f, :, i), view(g, :, i), advection, e[i]*Δt/Δp)
        end

        ey[k, :] = em.ey
        ex[k, :] = e
        ε_e[k] = integrate(x, e.^2)
        ε[k] = integrate(x, integrate(p, @. f*p^2)) + ε_e[k]
    end

    return (; t, x, p, n, ey, ex, ε_e, ε)
end
