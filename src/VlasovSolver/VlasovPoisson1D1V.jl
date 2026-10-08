"""
    VlasovPoisson1D1V

The 1D1V electrostatic driver: Strang-split Vlasov–Poisson on static grids,
with an optional collision operator, and the helpers it is built from.

It lived in `test/verification_harness.jl` until 0.2, which meant that every
number in the README's verification table was a statement about test code and
that no user could run a Landau damping problem without copying it.
"""
module VlasovPoisson1D1V

using ..Advection: AbstractAdvection1D, PFCNonUniform, OnGrid,
                   Upwind, Godunov, PiecewiseConstant, PiecewiseLinear, VanLeer, Superbee,
                   SemiLagrangian, LinearSpline, PFC
using ..StrangSplitting: strang_step!, Collide
using ..PoissonFourier1D: PoissonFFT1D, solve!
using LinearAlgebra: transpose!
import ..workspace

export vlasov_poisson, cell_widths, mode_amplitude, make_poisson

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
    p = uniform_poisson(x)
    ws = workspace(p)
    return (e, ρ) -> solve!(e, ρ, p, ws)
end

# `make_poisson`'s check and solver without the closure, for the driver, which
# keeps the solver and its workspace in its `FieldSolve`.
function uniform_poisson(x)
    length(x) ≥ 2 || throw(ArgumentError("the x grid needs at least 2 points, got $(length(x))"))
    Δx = x[2] - x[1]
    tol = max(1e-10*abs(Δx), 8eps(float(maximum(abs, x))))    # rounding; see above
    for j in 2:length(x)-1
        abs(x[j+1] - x[j] - Δx) ≤ tol || throw(ArgumentError(
            "the Poisson solve needs a uniform x grid: spacing $(x[j+1] - x[j]) at node $j, $Δx at node 1"))
    end
    return PoissonFFT1D(length(x), Δx)
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
keeps_bounds(og::OnGrid, lo, hi) = keeps_bounds(og.scheme, lo, hi)
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
    each is put on its grid by [`OnGrid`](@ref), which takes the step as a
    displacement and splits one wider than a cell.
  * `nᵢ`: the ion density over `x`; by default the Maxwellian's `Σ M Δv`,
    uniform. The ions are used as given, in the cell widths' type, and never
    rescaled.
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

**The run is in the data's type, `T = float(eltype(f₀))`.** `f`, the
field and every history are `T`, `E_modes` is `Complex{T}`: Float32 `f₀` gives
Float32 histories whatever the type of `t`, and Float64 `f₀` Float64 ones. The
cell widths are `T`, or a grid's own type where that is wider, and the density
`Σ f Δv`, the ions and the charge are in the widths' type: Float32 data on
Float64 grids solves for its field from a Float64 charge. The lines alone are
advanced in Float64 at least, as the default schemes are built, and rounded
into `f` once a step.

A run carries nothing from one step to the next but `f`, so it can be continued:
called again from the `f` it returned, over the rest of `t`, with the same ions,
schemes that do not depend on the starting `f`, and `renormalize = false`, it
takes the same steps to the bit. (Without `collisions` the defaults are bounded
above by the starting `f`, and the renormalisation is a factor of 1 only to
round-off.)

A scheme given for only one of `scheme_x`, `scheme_v` has to keep `f` inside the
default's bounds, `[0, maximum(f)]` (`[0, Inf)` with `collisions`), as
[`keeps_bounds`](@ref) decides, or the
call throws an `ArgumentError` before the first step. Pass both schemes, or a
`PFCNonUniform` with bounds both sweeps keep (`fmin = -Inf, fmax = Inf` for
none), for the other direction. With `collisions` the defaults have no upper
bound, since a collision step can raise a line's peak above `maximum(f)`.

Why each of these choices was made, with the measurements behind them, is in
`docs/src/driver-notes.md`. The run is [`setup`](@ref), then [`step!`](@ref) and
[`record!`](@ref) for each step.
"""
function vlasov_poisson(x, v, f₀, t;
                        scheme_x = nothing, scheme_v = nothing, invariants = false,
                        modes = (), nᵢ = nothing, renormalize = nᵢ === nothing,
                        collisions = nothing)
    prob, state, hist = setup(x, v, f₀, t; scheme_x, scheme_v, invariants, modes, nᵢ,
                              renormalize, collisions)
    run!(hist, state, prob, t)
    (; ε_e, ε, mass, momentum, l2, entropy, fmin, fmax, E_modes) = hist
    # `f` comes back too. It costs nothing -- the array exists either way -- and
    # it is the only way to ask a question about the distribution rather than
    # about a moment of it, which is what the reversibility test needs.
    return (; ε_e, ε, mass, momentum, l2, entropy, fmin, fmax, E_modes, f = state.f)
end

"""
    Problem

What a run of [`vlasov_poisson`](@ref) is, fixed for the whole run: the grids
`x` and `v` as given, their cell widths `Δx` and `Δv`, the ion density `nᵢ`,
the two schemes on their grids, `ox` and `ov` ([`OnGrid`](@ref)s), the collision
operator or `nothing`, and the wavenumbers of the modes recorded. Built by
[`setup`](@ref).

`W` is the widths' type, the data's or wider. The ions are `W`, as the density
they are set against is.
"""
struct Problem{W<:AbstractFloat, X<:AbstractVector, V<:AbstractVector,
               SX<:OnGrid, SV<:OnGrid, C, M}
    x::X
    v::V
    Δx::Vector{W}
    Δv::Vector{W}
    nᵢ::Vector{W}
    ox::SX
    ov::SV
    collisions::C
    modes::M
end

"""
    FieldSolve

The field solve of a step, called by [`strang_step!`](@ref) as its `cv`: from
`ft = f[v, x]` it sums the density `nₖ = Σ f Δv` over each column, solves for
the field `e` with `∂e/∂x = nₖ − nᵢ`, and returns the kick's displacements
`e·Δt`, all into buffers it owns. Mutable for `Δt` alone, which
[`step!`](@ref) sets before each step.

The density, the ions and the charge are in the widths' type `W`, as
`ft .* Δv` is: Float32 data on Float64 grids sums its density and cancels it
against the ions in Float64. Summed and cancelled in Float32, the field of a
Float32 Landau run on Float64 grids was 13 times further from the Float64 run's
than it had been; see `docs/src/driver-notes.md`.

`nrow` is `nₖ` as a `1 × Nx` matrix, sharing its memory, for `sum!`: that sums
each column in the order `sum(ft .* Δv, dims = 1)` does, to the bit, where a
loop or a matrix-vector product rounds differently. `tmp` holds `ft .* Δv`:
the state's scratch, the size of `f`, which the field solve and the diagnostics
take in turn, when `W` is `T`, and a matrix of the field solve's own when `W`
is wider.
"""
mutable struct FieldSolve{T<:AbstractFloat, W<:AbstractFloat, P, PW}
    Δt::T
    const e::Vector{T}
    const nₖ::Vector{W}
    const nrow::Matrix{W}
    const ρ::Vector{W}
    const αv::Vector{T}
    const nᵢ::Vector{W}
    const Δv::Vector{W}
    const tmp::Matrix{W}
    const poisson::P
    const pws::PW
end

function (fs::FieldSolve)(ft)
    @. fs.tmp = ft*fs.Δv
    sum!(fs.nrow, fs.tmp)
    @. fs.ρ = fs.nₖ - fs.nᵢ
    solve!(fs.e, fs.ρ, fs.poisson, fs.pws)
    @. fs.αv = fs.e*fs.Δt
    return fs.αv
end

"""
    State

What a run of [`vlasov_poisson`](@ref) carries from step to step, and the
buffers it is stepped in: `f` as `f[v, x]`, the layout every sum is taken over,
and `fxv`, the same state as `f[x, v]`, which [`strang_step!`](@ref) steps; the x
displacements `αx`; the scratch `tmp`, shared with `field` when the widths are
the data's type; the quadrature weights `wt = Δv .* Δx'`; the splitting's
workspace `ws`; the [`FieldSolve`](@ref) `field`; and `kinetic`, the kinetic energy `Σ f v² ΔvΔx` of `f` as it stands,
which `ε` is centred with.
"""
struct State{T<:AbstractFloat, W<:AbstractFloat, WS, F<:FieldSolve}
    f::Matrix{T}
    fxv::Matrix{T}
    αx::Vector{T}
    tmp::Matrix{T}
    wt::Matrix{W}
    ws::WS
    field::F
    kinetic::Base.RefValue{T}
end

"""
    Histories

The histories [`vlasov_poisson`](@ref) returns, one entry per time: `ε_e` and
`ε` always; `mass`, `momentum`, `l2`, `entropy`, `fmin` and `fmax` with
`invariants`, `nothing` without; `E_modes`, `Nt × length(modes)`, with `modes`.
All in the data's type `T`, `E_modes` in `Complex{T}`. Written by
[`record!`](@ref).
"""
struct Histories{T<:AbstractFloat, I<:Union{Nothing, Vector{T}}, E<:Union{Nothing, Matrix{Complex{T}}}}
    ε_e::Vector{T}
    ε::Vector{T}
    mass::I
    momentum::I
    l2::I
    entropy::I
    # `f` is a distribution function: negative values are unphysical, and a
    # scheme that produces them can be accurate on a linear diagnostic while
    # being unusable on a nonlinear run. Tracked because the comparison study
    # would otherwise rank a non-positive scheme first without saying so.
    fmin::I
    # And the maximum, which Liouville's theorem forbids to grow just as it
    # forbids the minimum to fall; the defaults hold their limiter to the
    # initial one, and this is how a run shows whether it stayed there.
    fmax::I
    # Per-mode field amplitudes, when asked for. `ε_e` sums every mode in the
    # box, which is fine while one of them dominates and misleading the moment
    # another does -- the recurrence of the second harmonic arrives at half the
    # time the seeded mode's does, and in `ε_e` it is indistinguishable from the
    # seeded mode coming back early.
    E_modes::E
end

"""
    setup(x, v, f₀, t; scheme_x, scheme_v, invariants, modes, nᵢ, renormalize, collisions)

A run of [`vlasov_poisson`](@ref) with the same arguments, before its first
step: `(prob, state, hist)`, a [`Problem`](@ref), a [`State`](@ref) holding `f₀`
as the run starts from it (rescaled if `renormalize`), and [`Histories`](@ref)
for `length(t)` entries. Every argument is checked here, and every buffer the
steps use is allocated here.
"""
function setup(x, v, f₀, t; scheme_x = nothing, scheme_v = nothing, invariants = false,
               modes = (), nᵢ = nothing, renormalize = nᵢ === nothing, collisions = nothing)
    # Every history entry k is taken during step k, and the last one is copied
    # from the one before: a single time has no step to take one in.
    length(t) ≥ 2 || throw(ArgumentError(
        "t needs at least 2 times, the start and one step, got $(length(t))"))
    # The run is in the data's type. The widths follow it, as the collision
    # scratch does: a Float32 or a Rational grid under Float64 `f` takes neither
    # its sums nor its default schemes out of Float64.
    T = float(eltype(f₀))
    Δx₀, Δv₀ = cell_widths(x), cell_widths(v)
    W = promote_type(T, eltype(Δx₀), eltype(Δv₀))
    Δx, Δv = convert(Vector{W}, Δx₀), convert(Vector{W}, Δv₀)

    if nᵢ === nothing
        # By the same sum over `v` as the electron density the field is solved
        # from, so that a uniform `f` is neutral as it stands.
        ions = fill(W(sum(@. exp(-0.5*v^2)/sqrt(2π)*Δv)), length(x))
    else
        length(nᵢ) == length(x) || throw(DimensionMismatch(
            "nᵢ has $(length(nᵢ)) points, the x grid $(length(x))"))
        ions = convert(Vector{W}, float.(nᵢ))
    end
    Nᵢ = sum(ions .* Δx)

    f = Matrix{T}(f₀)
    renormalize && (f .*= Nᵢ/sum(f .* (Δv .* Δx')))

    # Collisions relax a line towards a Maxwellian whose peak can sit above the
    # line's own, so the upper bound is not the run's to keep: the defaults then
    # bound `f` below only.
    bound = collisions === nothing ? maximum(f) : Inf
    # The defaults are Float64 at least, as the step works the line (`work`
    # below): a Float32 one rounds each flux difference, and can leave a cell an
    # ulp outside its bounds, which the next checked call refuses -- one
    # subnormal below 0 at the edge of a plateau, as `docs/src/driver-notes.md`
    # measures.
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

    # Each scheme on its grid, taking a displacement and splitting a step wider
    # than a cell; a uniform scheme on a non-uniform grid is refused here.
    ox, ov = OnGrid(sx, Δx), OnGrid(sv, Δv)
    # The step runs on `f[x, v]`, whose columns are the x lines, and sweeps v over
    # the transpose it keeps. The lines are worked in Float64 at least and
    # rounded into `f` once a step, whatever `f`'s own type; see the defaults.
    fxv = permutedims(f)
    work = promote_type(Float64, T)
    ws = workspace(ox, ov, fxv, collisions, work)

    Nx, Nv = length(x), length(v)
    poisson = uniform_poisson(x)
    # One scratch matrix, reused: the density, the energy and the four
    # invariants are sums over arrays the shape of `f`, and allocating a
    # temporary per sum per step dominated the step itself.
    tmp = similar(f)
    nₖ = Vector{W}(undef, Nx)
    field = FieldSolve(zero(T), Vector{T}(undef, Nx), nₖ, reshape(nₖ, 1, Nx),
                       Vector{W}(undef, Nx), Vector{T}(undef, Nx), ions, Δv,
                       W === T ? tmp : similar(f, W), poisson, workspace(poisson))
    wt = Δv .* Δx'      # the flux-form quadrature, for every sum
    @. tmp = f*(v*v)*wt     # as `record!` takes it
    state = State(f, fxv, Vector{T}(undef, Nv), tmp, wt, ws, field, Ref(sum(tmp)))

    Nt = length(t)
    history() = Vector{T}(undef, Nt)
    invariant() = invariants ? history() : nothing
    hist = Histories(history(), history(), invariant(), invariant(), invariant(),
                     invariant(), invariant(), invariant(),
                     isempty(modes) ? nothing : zeros(Complex{T}, Nt, length(modes)))
    return Problem(x, v, Δx, Δv, ions, ox, ov, collisions, modes), state, hist
end

"""
    step!(state, prob, Δt)

One Strang-split step of `state.f` by `Δt`: [`strang_step!`](@ref) on
`state.fxv`, with the field solved by `state.field` after the first half step
and the collisions, if any, as its [`Collide`](@ref) hook; then `f[v, x]` from
the result. Allocates nothing that grows with the grid.
"""
function step!(s::State, p::Problem, Δt)
    s.field.Δt = Δt
    @. s.αx = p.v*Δt
    hook = p.collisions === nothing ? nothing : Collide(p.collisions, p.v, Δt)
    strang_step!(s.fxv, p.ox, p.ov, s.αx, s.field, s.ws, hook)
    # Back to `f[v, x]`, so that every sum below runs in the order it always
    # has: summed over `[x, v]` the same terms round differently.
    transpose!(s.f, s.fxv)
    return s
end

"""
    record!(hist, k, state, prob)

Entry `k` of each history from the step just taken: the field's energy and
modes as the step solved them, `ε` with the kinetic energy centred on the kick,
and the invariants of `f` after the step. Allocates nothing.
"""
function record!(h::Histories{T}, k, s::State, p::Problem) where {T}
    e, f, tmp, wt, v, Δx = s.field.e, s.f, s.tmp, s.wt, p.v, p.Δx
    if h.E_modes !== nothing
        for (j, km) in enumerate(p.modes)
            h.E_modes[k, j] = mode_amplitude(e, p.x, T(km))
        end
    end
    h.ε_e[k] = sum(e[j]^2*Δx[j] for j in eachindex(e))
    # The kick is all that moves the kinetic energy: the x half-steps on
    # either side of it carry each row of `f` in flux form, which keeps the
    # row's sum and so `Σ f v² ΔvΔx`. Its value when `e` was solved is then
    # the mean of the two sides of the kick, and `ε` is a sample of the same
    # instant as `ε_e`; see `vlasov_poisson`. `v*v` is `v^2` to the bit, but
    # the literal power boxed a value every step on Julia 1.10 under
    # `--check-bounds=yes`, which `Pkg.test` passes.
    @. tmp = f*(v*v)*wt
    before = s.kinetic[]
    s.kinetic[] = sum(tmp)
    h.ε[k] = (before + s.kinetic[])/2 + h.ε_e[k]

    if h.mass !== nothing
        @. tmp = f*wt
        h.mass[k] = sum(tmp)
        @. tmp = f*v*wt
        h.momentum[k] = sum(tmp)
        @. tmp = f*f*wt
        h.l2[k] = sum(tmp)
        @. tmp = nlogn(f)*wt
        h.entropy[k] = sum(tmp)
        h.fmin[k] = minimum(f)
        h.fmax[k] = maximum(f)
    end
    return h
end

# The time loop, behind a function barrier: what `setup` returns has types that
# depend on the keywords' values, and this is compiled for the ones it got.
function run!(h::Histories, s::State, p::Problem, t)
    for k in 1:length(t)-1
        step!(s, p, t[k+1] - t[k])
        record!(h, k, s, p)
    end
    for y in (h.ε_e, h.ε, h.mass, h.momentum, h.l2, h.entropy, h.fmin, h.fmax)
        y === nothing || (y[end] = y[end-1])
    end
    h.E_modes === nothing || (@views h.E_modes[end, :] .= h.E_modes[end-1, :])
    return h
end

end # module
