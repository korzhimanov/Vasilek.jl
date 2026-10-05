"""
    VlasovPoisson1D1V

The 1D1V electrostatic driver: Strang-split Vlasov–Poisson on static grids,
with an optional collision operator, and the helpers it is built from.

It lived in `test/verification_harness.jl` until 0.2, which meant that every
number in the README's verification table was a statement about test code and
that no user could run a Landau damping problem without copying it.
"""
module VlasovPoisson1D1V

using ..Advection: AbstractAdvection1D, PFCNonUniform, advect!,
                   Upwind, Godunov, PiecewiseConstant, PiecewiseLinear, VanLeer, Superbee,
                   SemiLagrangian, LinearSpline, PFC
using ..StrangSplitting: make_time_step_2d!
using ..Collisions: collide!
using ..PoissonFourier1D: PoissonFFT1D, solve!
import ..workspace

export vlasov_poisson, cell_widths, line_advector, substeps, mode_amplitude, make_poisson

"""
    make_poisson(x)

Spectral Poisson solve on a uniform periodic `x` grid, returning `(e, ρ) -> e`
with `∂e/∂x = ρ`: [`PoissonFFT1D`](@ref) with its centred derivative and a
workspace of its own.

The grid is checked to be uniform to rounding: each spacing must match the
first to 1e-10 of it or to 8 ulps of the largest coordinate, whichever is more.
A stretched grid would otherwise get a field from its first spacing alone.
"""
function make_poisson(x)
    length(x) ≥ 2 || throw(ArgumentError("the x grid needs at least 2 points, got $(length(x))"))
    Δx = x[2] - x[1]
    tol = max(1e-10*abs(Δx), 8eps(float(maximum(abs, x))))    # rounding; see above
    for j in 2:length(x)-1
        abs(x[j+1] - x[j] - Δx) ≤ tol || throw(ArgumentError(
            "the Poisson solve needs a uniform x grid: spacing $(x[j+1] - x[j]) at node $j, $Δx at node 1"))
    end
    p = PoissonFFT1D(length(x), Δx)
    ws = workspace(p)
    return (e, ρ) -> solve!(e, ρ, p, ws)
end

"""
    cell_widths(z)

The width each node of the grid `z` stands for: the spacing at the two ends,
and half the distance between its neighbours inside. These are the weights of
the sums `Σ f ΔvΔx` that the flux-form schemes conserve. A floating-point grid
gets widths of its own element type.
"""
cell_widths(z) = vcat([z[2]-z[1]], (z[3:end] - z[1:end-2])/2, [z[end]-z[end-1]])


"""
    nlogn(u)

`-u·log u`, and zero at `u ≤ 0`, for the entropy integrand. A scalar function
rather than an `ifelse` in a broadcast, which would evaluate `log` of the
negative values schemes without a positivity limiter produce.
"""
nlogn(u) = u > 0 ? -u*log(u) : zero(u)

"""
    line_advector(scheme, Δz)

`scheme` as an in-place `(column, α)` advector, where `α` is always a
**displacement**, whatever the scheme's own fourth argument is. Uniform-grid
schemes get `α/Δz`, and are refused on a non-uniform grid rather than given one
of its spacings. A displacement wider than the narrowest cell is split into the
fewest equal sub-steps that fit ([`substeps`](@ref)) for `PFCNonUniform`; a
uniform scheme's `advect!` refuses it.

The line is worked in Float64, its buffer and the scheme's scratch alike, and
rounded into the column once per step, whatever the scheme's own element type.
"""
function line_advector(scheme::PFCNonUniform, Δz)
    n = length(Δz)
    # Not the scheme's type: a Float32 scheme's accumulator rounded a Float64
    # line to Float32 on every sub-step.
    ws = workspace(scheme, n, Float64)
    buf = Vector{Float64}(undef, n)
    h = minimum(Δz)
    return function (col, α)
        m = substeps(α, h)
        for _ = 1:m
            advect!(buf, col, scheme, α/m, ws)
            copyto!(col, buf)
        end
        return col
    end
end

"""
    substeps(α, h)

The fewest equal parts of `α` that are each no longer than `h`, and one for an
`α` that is not finite, which `advect!` then refuses itself.
"""
function substeps(α, h)
    isfinite(α) || return 1
    m = max(1, ceil(Int, abs(α)/h))
    return abs(α/m) > h ? m + 1 : m    # the quotient can round below the integer
end

function line_advector(scheme, Δz)
    n = length(Δz)
    h = Δz[1]
    all(d -> isapprox(d, h; rtol = 1e-12), Δz) || error(
        "$(nameof(typeof(scheme))) takes a Courant number, which a non-uniform grid " *
        "does not have (spacings range over $(extrema(Δz))); use PFCNonUniform here")
    ws = workspace(scheme, n, Float64)
    buf = Vector{Float64}(undef, n)
    return (col, α) -> (advect!(buf, col, scheme, α/h, ws); copyto!(col, buf))
end

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
    keeps_bounds(scheme, lo, hi)

Whether `scheme` keeps data that starts inside `[lo, hi]` inside it, to
round-off: what [`vlasov_poisson`](@ref) asks of a scheme given for one
direction while the other is left to its default, a `PFCNonUniform` built on
`[0, maximum(f)]` that refuses a line outside them.

True for the schemes that keep their data between its own extrema -- `Upwind`,
`Godunov` with a constant reconstruction or a limiter, the linear
`SemiLagrangian` -- and for a `PFC` or `PFCNonUniform` whose own bounds lie
inside `[lo, hi]`. False for `LaxWendroff`, `Godunov(PiecewiseLinear())` without
a limiter and the quadratic and cubic `SemiLagrangian`, which are linear and
above first order and so, by Godunov's theorem, overshoot; for a `PFC` bounded
wider, whose limiter lets the data out to its own bounds; and for a scheme of any
other type, which the driver cannot vouch for. The measurements behind it are in
`docs/src/driver-notes.md`.
"""
keeps_bounds(::AbstractAdvection1D, lo, hi) = false
keeps_bounds(::Upwind, lo, hi) = true
keeps_bounds(::Godunov{PiecewiseConstant}, lo, hi) = true
keeps_bounds(::Godunov{PiecewiseLinear, <:Union{VanLeer, Superbee}}, lo, hi) = true
keeps_bounds(::SemiLagrangian{LinearSpline}, lo, hi) = true
# In the scheme's own precision, `[lo, hi]` rounded outward as its constructor
# rounds its bounds: a Float32 scheme built on `maximum(f)` keeps it.
function keeps_bounds(s::Union{PFC, PFCNonUniform}, lo, hi)
    T = typeof(s.fmax)
    return T(lo, RoundDown) ≤ s.fmin && s.fmax ≤ T(hi, RoundUp)
end

# A scheme given for one direction while the other is left to its default has to
# keep `f` inside the default's bounds; see `vlasov_poisson`.
function refuse_unbounded(given, name, other, bound)
    keeps_bounds(given, 0.0, bound) && return nothing
    throw(ArgumentError(
        "$name = $(sprint(show, given)) may take f outside [0, $bound], " *
        "the bounds of the PFCNonUniform that $other defaults to, which stops the run at " *
        "the first line handed to it outside them. Pass $other as well: the same scheme, " *
        "or a PFCNonUniform with bounds both sweeps keep (fmin = -Inf, fmax = Inf for none)"))
end

"""
    vlasov_poisson(x, v, f₀, t; scheme_x, scheme_v, invariants = false, modes = (),
                   nᵢ = nothing, renormalize = nᵢ === nothing, collisions = nothing)

Strang-split electrostatic Vlasov–Poisson for electrons over fixed ions, from
`f₀` (indexed `f[v, x]`) through the times `t`; Vlasov–Poisson–BGK when
`collisions` is a collision operator. Units are those of
`docs/src/normalization.md`.

  * `x`: a periodic grid, uniform (checked). `v`: any grid, uniform or not.
  * `t`: at least two times; steps may differ.
  * `scheme_x`, `scheme_v`: advection schemes, or functions of the starting `f`
    that return one. Both default to `PFCNonUniform` on their grid, in Float64
    at least, bounded by `[0, maximum(f)]` (below only with `collisions`);
    their steps wider than a cell are split by [`line_advector`](@ref).
  * `nᵢ`: the ion density over `x`; by default the Maxwellian's `Σ M Δv`,
    uniform. The ions are used as given, and never rescaled.
  * `renormalize`: rescale `f` on entry by a single factor, so that its charge
    `Σ f ΔvΔx` equals the ions' `Σ nᵢ Δx`. On by default without `nᵢ` and off
    with it. With a given `nᵢ`, `renormalize = true` matches the totals only:
    the box becomes neutral as a whole, and `f` keeps its shape over `x`
    rather than being fitted to the ions' profile.
  * `collisions`: applied to every velocity line for half a step either side of
    the kick, `X(Δt/2)·C(Δt/2) K(Δt) C(Δt/2)·X(Δt/2)`, so the step stays second
    order.
  * `modes = (k₁, k₂, …)`: record each mode's complex field amplitude,
    [`mode_amplitude`](@ref).

Returns a NamedTuple:

  * `ε_e`, `ε`: electric and total energy, the cell-width sums `Σ E²Δx` and
    `Σ f v² ΔvΔx + Σ E²Δx` (twice the physical energies);
  * `E_modes` (with `modes`), `Nt × length(modes)`;
  * `mass`, `momentum`, `l2`, `entropy`, `fmin`, `fmax` (with `invariants`),
    as cell-width sums;
  * `f`, the final distribution.

**Row `k` of `E_modes`, `ε_e` and `ε` is sampled at `t[k] + Δt/2`, not `t[k]`;
row `k` of `mass`, `momentum`, `l2`, `entropy`, `fmin` and `fmax` at the end of
step `k`, `t[k+1]`.** None is sampled at `t[1]`, and the last row of each is
a copy of the one before. The kinetic energy in `ε` is centred on the kick, so
`ε` and `ε_e` describe the same instant.

A scheme given for only one of `scheme_x`, `scheme_v` has to keep `f` inside the
default's bounds, `[0, maximum(f)]` (`[0, Inf)` with `collisions`), as
[`keeps_bounds`](@ref) decides, or the
call throws an `ArgumentError` before the first step. Pass both schemes, or a
`PFCNonUniform` with bounds both sweeps keep (`fmin = -Inf, fmax = Inf` for
none), for the other direction. With `collisions` the defaults have no upper
bound, since a collision step can raise a line's peak above `maximum(f)`.

Why each of these choices was made, with the measurements behind them, is in
`docs/src/driver-notes.md`.
"""
function vlasov_poisson(x, v, f₀, t;
                        scheme_x = nothing, scheme_v = nothing, invariants = false,
                        modes = (), nᵢ = nothing, renormalize = nᵢ === nothing,
                        collisions = nothing)
    # Every history entry k is taken during step k, and the last one is copied
    # from the one before: a single time has no step to take one in.
    length(t) ≥ 2 || throw(ArgumentError(
        "t needs at least 2 times, the start and one step, got $(length(t))"))
    # The widths follow the data, as the collision scratch does: a Float32 or a
    # Rational grid under Float64 `f` takes neither its sums nor its default
    # schemes out of Float64.
    Δx₀, Δv₀ = cell_widths(x), cell_widths(v)
    W = promote_type(float(eltype(f₀)), eltype(Δx₀), eltype(Δv₀))
    Δx, Δv = convert(Vector{W}, Δx₀), convert(Vector{W}, Δv₀)

    if nᵢ === nothing
        # By the same sum over `v` as the electron density the field is solved
        # from, so that a uniform `f` is neutral as it stands.
        nᵢ = fill(sum(@. exp(-0.5*v^2)/sqrt(2π)*Δv), length(x))
    else
        length(nᵢ) == length(x) || throw(DimensionMismatch(
            "nᵢ has $(length(nᵢ)) points, the x grid $(length(x))"))
        nᵢ = collect(float.(nᵢ))
    end
    Nᵢ = sum(nᵢ .* Δx)

    f = copy(f₀)
    renormalize && (f .*= Nᵢ/sum(f .* (Δv .* Δx')))
    g = f'

    # Collisions relax a line towards a Maxwellian whose peak can sit above the
    # line's own, so the upper bound is not the run's to keep: the defaults then
    # bound `f` below only.
    bound = collisions === nothing ? maximum(f) : Inf
    # The defaults are Float64 at least, as `line_advector` works the line: a
    # Float32 one rounds each flux sum, and at a line's peak the sum can land an
    # ulp above `fmax`, which the next checked call refuses.
    default(widths) = PFCNonUniform(convert(Vector{promote_type(Float64, eltype(widths))}, widths);
                                    fmin = 0.0, fmax = bound)
    pick(s, widths) = s === nothing ? default(widths) :
                      s isa AbstractAdvection1D ? s : s(f)
    sx, sv = pick(scheme_x, Δx), pick(scheme_v, Δv)
    # A default stops the run at the first line it is handed outside its bounds;
    # a scheme given for the other direction that need not keep them is refused
    # here instead, before the first step. See the docstring.
    scheme_x === nothing && scheme_v !== nothing &&
        refuse_unbounded(sv, :scheme_v, :scheme_x, bound)
    scheme_v === nothing && scheme_x !== nothing &&
        refuse_unbounded(sx, :scheme_x, :scheme_v, bound)

    advect_x! = line_advector(sx, Δx)
    advect_v! = line_advector(sv, Δv)
    if collisions !== nothing
        cws = workspace(collisions, length(v), eltype(f))
        cbuf = similar(v)
    end

    solve_poisson! = make_poisson(x)
    e = similar(x)
    ε_e = similar(t)
    ε = similar(t)
    # One scratch matrix, reused: the energy and the four invariants are five
    # integrals of the same shape as `f`, and allocating a temporary per integral
    # per step dominated the step itself when this was first written.
    tmp = similar(f)
    wt       = Δv .* Δx'      # the flux-form quadrature, for every sum below
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
    @. tmp = f*v^2*wt
    kinetic = sum(tmp)        # before the first kick

    for k in 1:length(t)-1
        Δt = t[k+1] - t[k]
        vΔt(_) = v*Δt
        function eΔt(ff)
            nₖ = vec(sum(ff'.*Δv, dims = 1))
            solve_poisson!(e, nₖ - nᵢ)
            return e*Δt
        end
        kick! = collisions === nothing ? advect_v! : function (col, α)
            collide!(cbuf, col, collisions, v, Δt/2, cws); copyto!(col, cbuf)
            advect_v!(col, α)
            collide!(cbuf, col, collisions, v, Δt/2, cws); copyto!(col, cbuf)
            return col
        end
        make_time_step_2d!((g, f), (vΔt, eΔt), (advect_x!, kick!))
        if E_modes !== nothing
            for (j, km) in enumerate(modes)
                E_modes[k, j] = mode_amplitude(e, x, km)
            end
        end
        ε_e[k] = sum(e[j]^2*Δx[j] for j in eachindex(e))
        # The kick is all that moves the kinetic energy: the x half-steps on
        # either side of it carry each row of `f` in flux form, which keeps the
        # row's sum and so `Σ f v² ΔvΔx`. Its value when `e` was solved is then
        # the mean of the two sides of the kick, and `ε` is a sample of the same
        # instant as `ε_e`; see the docstring.
        @. tmp = f*v^2*wt
        before, kinetic = kinetic, sum(tmp)
        ε[k] = (before + kinetic)/2 + ε_e[k]

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
    end
    for h in (ε_e, ε, mass, momentum, l2, entropy, fmin, fmax)
        h === nothing || (h[end] = h[end-1])
    end
    E_modes === nothing || (E_modes[end, :] = E_modes[end-1, :])
    # `f` comes back too. It costs nothing -- the array exists either way -- and
    # it is the only way to ask a question about the distribution rather than
    # about a moment of it, which is what the reversibility test needs.
    return (; ε_e, ε, mass, momentum, l2, entropy, fmin, fmax, E_modes, f)
end

end # module
