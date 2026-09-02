# Shared 1D1V Vlasov–Poisson driver for the extended verification tests.
# Deliberately close to what the verification notebooks do, so the two agree.

using Vasilek
using Vasilek: StrangSplitting
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
-1.3e-10 and -4.4e-10 on the large-amplitude case in
`verification/scheme-comparison.jl`, and the entropy diagnostic threw
`DomainError` on both until this was split out.
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
