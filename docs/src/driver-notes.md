# Notes on the 1D1V driver

```@meta
CurrentModule = Vasilek.VlasovPoisson1D1V
```

The reasoning behind `Vasilek.VlasovPoisson1D1V`, moved here from its
docstrings so that those can say what the functions do. The numbers below are
as they were measured when each choice was made; the test suite asserts the
ones that matter.

## `make_poisson`

Spectral Poisson solve on a uniform periodic x grid, returning `e` with
`∂e/∂x = ρ`: the package's `PoissonFFT1D` with its centred derivative and a
workspace of its own.

The grid must be uniform, and is checked: the solve takes its spacing from the
first two points, so a stretched x grid gave a wrong field with no error (16% in
the field energy for spacings within ±8% of the mean). The advection would take
one; the field solve would need a different method.

Uniform is to rounding: each spacing must match the first to 1e-10 of it, or
to 8 ulps of the largest coordinate where that is more. The points are rounded
on the scale of the largest of them, so the spacings of a uniform grid scatter
by an ulp or two of it whatever built it -- at most 2 for a collected range,
`LinRange`, `(1:N)*Δx` or `cumsum`, measured to N = 10⁶ in Float64 and in
Float32. Against the spacing that is N·eps: far below 1e-10 for a Float64 grid
of ordinary size, but a few parts in 10⁶ for a 64-point Float32 one, and past
1e-10 for Float64 ones of N ≈ 5×10⁵ or far from the origin. The ulps admit
those; a grid stretched by a millionth of a cell is refused either way.

This used to be a second implementation, which took its wavenumbers from a
float range whose length rounds: on the 8π box 61 of the even `Nx` in 8:512
came out one short and the solve threw on its first call. The verification
runs now exercise the solver the package ships.

## `nlogn`

    nlogn(u)

`-u·log u`, and zero at `u ≤ 0`, for the entropy integrand.

A scalar function rather than an `ifelse` inside the broadcast, because
`ifelse` is an ordinary call and evaluates **both** arguments: written that way
the guard does not guard, and `log` is handed the negative value anyway. That is
not hypothetical here -- `LaxWendroff` and cubic `SemiLagrangian` drive `f` to
-0.093 and -0.105 on the large-amplitude case in
`verification/scheme-comparison.jl`, against a peak of 0.6, and the entropy
diagnostic threw `DomainError` on both until this was split out. That is an 18%
undershoot of the peak rather than round-off leaking below zero, which is worth
stating precisely: it is the size of the overshoot that makes the guard a
statement about the schemes rather than about floating point.

## `line_advector`

    line_advector(scheme, Δz)

Wrap `scheme` as an in-place `(column, α)` advector, where `α` is always a
**displacement** -- a length -- whatever the scheme's own fourth argument means.

`PFCNonUniform` takes a displacement already. Every other scheme takes a Courant
number, which only exists on a uniform grid, so the wrapper divides by the
spacing for those and **refuses** a non-uniform grid rather than picking one of
its spacings and being quietly wrong by the ratio between them. That asymmetry
is a documented wart of the advection API (see
[normalization.md](normalization.md)); this is the one place the verification
runs have to absorb it.

**A displacement wider than the narrowest cell is split** into the fewest equal
sub-steps that fit ([`substeps`](@ref)). `advect!` refuses it whole, and a
translation by `α` is `m` translations by `α/m`. Every call the suite makes
through here is within the bound, and so one call exactly as before, except in
the two-stream runs, where the field drives the velocity sweep: they keep
growing after their fits, into saturation, and the fastest would ask for 217
cells a step by `tmax`. That moves no measured number, since no fit contains a
split step. The strong-damping run on the notebook's non-uniform grid comes
closest otherwise: its field peaks at 0.9938 on the first step, with `Δt` equal
to the narrow cells' width. It crossed, at 1.0017, while the driver still
rescaled `f` by the trapezoid, 0.79% up at that amplitude.

The uniform-grid method does not split. Nothing here asks a uniform scheme for
more than `c = 0.50`, and `advect!` says so if something ever does.

## `vlasov_poisson`

    vlasov_poisson(x, v, f₀, t; scheme_x, scheme_v, invariants = false, modes = (),
                   nᵢ = nothing, renormalize = nᵢ === nothing, collisions = nothing)

Strang-split Vlasov–Poisson on a static grid, and Vlasov--Poisson--BGK when
`collisions` is a collision operator.

Returns a NamedTuple with `ε_e` (electric energy) and `ε` (total energy)
histories. With `invariants = true` it also returns the `mass`, `momentum`,
`l2` and `entropy` histories -- the conserved quantities that are *not* the
energy, and that nothing asserted until now. `modes = (k₁, k₂, …)` adds
`E_modes`, an `Nt × length(modes)` matrix of complex field amplitudes from
[`mode_amplitude`](@ref), which is what separates a mode from its harmonics
where `ε_e` cannot.

**Row `k` of `E_modes`, `ε_e` and `ε` is sampled at `t[k] + Δt/2`, not `t[k]`;
row `k` of `mass`, `momentum`, `l2`, `entropy`, `fmin` and `fmax` at the end of
step `k`, `t[k+1]`.** None is sampled at `t[1]`, and the last row of each is
a copy of the one before.
The field is solved inside the step, after Strang's first `x` half-step, and
recorded from there. A rate or a frequency cannot see a constant time offset;
a *phase* can, and does: removing a Doppler factor `exp(-ikut)` at `t[k]` leaves
`k·u·Δt/2` behind on every sample, which the Galilean test in
`test_verification.jl` once reported as the grid's own non-invariance.

**The invariants and the energies use the cell-width sums `Σ f ΔvΔx` and
`Σ e² Δx`, not `integrate`.** That is the quadrature the schemes actually
conserve: `PFC` is a flux form, so what leaves one cell enters its neighbour and
the full-weight sum is preserved exactly. The trapezoid halves the two endpoint
weights, which no flux conservation law protects, and measuring with it reports a
drift that belongs to the quadrature rather than to the scheme. Measured over the
k = 0.5 Landau run, 875 steps: mass drifts 1.6e-15 by the cell-width sum against
**2.4e-4** by the trapezoid, and momentum stays at 2.2e-15 against 7.3e-4. Both
trapezoid figures are the endpoint weighting, not the solver.

The energies kept the trapezoid longer, on the grounds that the tolerances they
meet are half a percent. For a standing wave that does little harm: the field
`∝ sin kx` has a node on the seam, `x = L ≡ 0`, where the trapezoid's two
half-weighted end points sit. A travelling wave crosses the seam, and the
trapezoid books its passage as energy. On the bump-on-tail instability of Arber
and Vann, seeded at 1e-6 on 64 × 361 cells, where `ε` is 61.8, it swung by
−0.346 to +0.175 over t ∈ [60, 100], and by −0.185 to +0.093 at Nx = 128 --
first order in Δx, where the sums here drift by 9.7e-3 and 1.0e-3. The x
half-steps alone, which cannot change the kinetic energy, moved the trapezoid's
by up to 3.9e-4 of it per half-step. The Galilean test's boosted mode travels
too, and its fitted damping rate read 0.136% off the rest frame's through the
same effect, where it now reads 0.014%.

**The kinetic energy in `ε` is centred on the kick.** It is summed after the
step, the field was solved before the kick, and the x half-steps in between
leave `Σ f v² ΔvΔx` alone -- by 2.4e-16 of it per half-step at worst, over the
bump-on-tail run -- so the kick is the only change, and the kinetic energy at
`t[k] + Δt/2` is the mean of its two sides. Summed after the kick, as it was, `ε`
counted half the kick's work early: an error first order in Δt that oscillates
with the power the field exchanges with the particles. On the plasma oscillation
the energy test runs, that alone swung `ε` by 1.2e-3 of itself within each plasma
period, and by 1.7e-3 with the trapezoid, against the 3.8e-3 it drifts through
t = 3000; centred, it swings by 6.6e-5.

**The ions are a fixed background**, a Maxwellian's density on the grid unless
`nᵢ` gives a profile over `x`. Without `nᵢ`, `f` is rescaled on entry so that the
two integrate to the same charge; with it, `f` is taken as it comes; and
`renormalize` overrides either. A caller that passes `nᵢ` is handing over the
ions it means, matched to `f` or deliberately not, and the driver keeps them.

**The charges are the cell-width sums too**, `Σ f ΔvΔx` and `Σ nᵢ Δx`, and the
default `nᵢ` is the Maxwellian's `Σ M Δv`, the same sum over `v` the field is
solved from. They were the trapezoid, which on the periodic grid weights the two
end points by half. That is exact for proportional profiles only, and
`M(v)(1 + α cos kx)` is not one: for it the trapezoid's rescaling was
`≈ 1 + α/(Nx − 1)`, 7.9e-3 at α = 0.5 on 64 cells, and `ωₚ²` with it, where the
sums give 1 to round-off. It also rescaled again on every restart, where the
flux form has kept `Σ f ΔvΔx` and the sums find nothing to do. On the matched
pair of the verification harness's `bgk_equilibrium` it was 1 − 1.9e-3 and
doubled the equilibrium's drift; it is now 1 to 3.8e-15 there as well.

`scheme_x` and `scheme_v` default to `PFCNonUniform` on the two grids, which is
what the verification notebooks use and what every previous caller got. They are
arguments so that the same driver can measure what the physics costs under a
*different* scheme, which is what `verification/scheme-comparison.jl` does, and
so that a refinement study can hold the scheme fixed while moving the grid.

**Either may also be a function of the initial `f` that returns a scheme**, as
in `f -> PFC(fmin = 0.0, fmax = maximum(f))`, and is called with the `f` the run
actually starts from. That is the only way for a caller to bound a scheme by
that distribution: the driver rescales what it is handed before it runs, to the
ions' charge, so a bound taken from `f₀` beforehand sits below the maximum of
any `f` it scales up and trips `PFC`'s own check on the first call.

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
`test_verification.jl` converges at second order (×4.3, then ×4.1) where it
converged at third (×5.7, ×7.1). Elsewhere the numbers moved little and are
updated where they are quoted. No run leaves its bound: the `fmax` history of a
Landau, a strong Landau, a two-stream and an equilibrium run tops out exactly at
it, or below, to the last bit. `PFCNonUniform` now checks every call it takes,
sub-steps included, and none of them is out of bounds -- the two-stream runs
carried past their velocity Courant limits into saturation among them, since
[`line_advector`](@ref) splits every step that would cross a cell.

**`collisions` puts a collision operator in the step**, applied to every velocity
line for half a step on either side of the kick:

    X(Δt/2) · C(Δt/2) K(Δt) C(Δt/2) · X(Δt/2)

The v-sweep is where a column of `f` is a velocity line at fixed `x`, which is
all a collision operator acts on, and the composition stays symmetric, so the
step stays second order. `BGK` conserves each line's density, momentum and
energy, so the field solved before the kick is still the field after the first
half-collision, and the kinetic energy is still moved by the kick alone -- the
centring above holds as it stands.

**A scheme given for one direction while the other is left to its default has
to keep `f` inside those bounds, or the call is refused**, with an
`ArgumentError`, before the first step. Without that check the default refused
such a call later, mid-run: it stops the run at the first line handed to it
outside its bounds, and a scheme that does not keep them hands it one soon
enough. At 50% amplitude on one wavelength of `k = 0.5`, 64 × 121 cells over ±6
(the case in `test_driver.jl`), `LaxWendroff` in `v` took `f` to −1.9e-9 on the
16th step and the cubic `SemiLagrangian` to −7.0e-10 on the 21st.
[`keeps_bounds`](@ref) decides. It passes the schemes that keep their data
between its own extrema -- `Upwind`, `Godunov` with a constant reconstruction or
a limiter, the linear `SemiLagrangian` -- and a `PFC` or `PFCNonUniform` bounded
inside `[0, maximum(f)]`, as `f -> PFC(fmin = 0.0, fmax = maximum(f))` is. It
refuses `LaxWendroff`, `Godunov(PiecewiseLinear())` without a limiter, the
quadratic and cubic `SemiLagrangian`, a `PFC` bounded wider, and a scheme of any
other type. Given both schemes, the driver checks neither: the same scheme for
both, as `verification/scheme-comparison.jl` passes them, or for the other
direction a `PFCNonUniform` with bounds both sweeps keep -- `fmin = -Inf,
fmax = Inf` for none, which is PFC's reconstruction with nothing to limit it.

*Superseded for the upper bound:* with `collisions` the defaults no longer carry
`maximum(f)` as an upper bound, so the `DomainError`s described next no longer
occur from the defaults; they are kept as the reason for dropping it.

The check was about the advection alone. `collisions` can take `f` past
`maximum(f)` by themselves, since BGK relaxes a line towards a Maxwellian whose
peak may sit above the line's own, and the default's upper bound then stops the
run mid-way whatever scheme was given, or none. With `Godunov(PiecewiseLinear(),
Superbee())` in `v`, which passes the check, and `BGK(1.0)` on the 50% case
above, the default stops the run on the 2nd step, `f` at 5.4e-8 of `maximum(f)`
above it (a run that nothing stopped would climb to 3.3e-7). With both defaults
and `BGK(0.1)`, a line flat over |v| < 2 is stopped on the 1st step at 0.4165
against a bound of 0.3846, 8.3% above. Both are a `DomainError` from inside the
step, as before this check.

Lifting the bound moved the collisional numbers slightly, though it had ended
none of those runs: tight, it clipped the reconstruction in the peak cell. At
ν = 0.1, 0.3 and 1, γ against the root on the grid's field went from +0.37%,
+0.48%, +0.75% to +0.37%, +0.49%, +0.76% (γ by 3e-5 at most); the heat mode
from 0.43383 to 0.43377 against the root's 0.43369; `PartialBGK{(:n, :u)}`'s ω
from +0.04% to −0.10% of its own root. The energy holds better: at ν = 1 it
drifts 8.0e-9 on ±8, where it drifted 1.4e-7, and 8.0e-8 on ±4, where it
drifted 2.1e-7. The sampled Maxwellian on ±4, which the bound stopped on the
first step and `verification/collisional-damping.jl` ran with 50% of headroom,
now runs under the defaults and cools by the same −5.3%.

That refuses calls that ran. At 1% on ±4, where `f` comes no nearer 0 than
1.3e-4, `LaxWendroff` and the cubic spline stay inside the bounds in either
direction, to the end of a run. That was the problem's doing rather than the
call's -- at 50% the same calls in `v` stop on steps 16 and 21 -- and a check
made before the run sees only the call. The two ways of letting them all run
were measured and are worse. A default that does not check its bounds keeps a
limiter that data outside them turns inside out: with `LaxWendroff` in `v`, `f`
ends 16% of the peak away from where the same run ends with no bounds at all,
at a minimum of −0.041 against −0.104 -- corrupted, and looking the better for
it. A default that drops its bounds runs, but it is not the scheme documented
here, and it changes what a call that ran measures: `scheme_x = LaxWendroff()`
at 1% reads γ 0.156% above the root with the default `v` bounded, and 0.289%
with its bounds dropped.

### `keeps_bounds`

Measured over 10080 runs of 20 steps on rough data in `[0, 1]` -- uniform
random, a quarter of it zeroed, a step down to 1e-9, spikes on a floor of
1e-300, a Maxwellian -- at twelve Courant numbers across `[-1, 1]`, every scheme
on the same data. Those it accepts never took the data below 0, and its maximum
moved by round-off alone: 7.3e-16 of it at worst for `Godunov`, in 7.5% of the
runs, and 6.7e-16 for `PFC` bounded at the data's own maximum, which is the
default's case. Those it refuses took the data below 0 in 62% to 64% of the
runs, by up to 0.33 of its maximum, and above the maximum in 29% to 45%, by up
to 0.24; a `PFC` bounded at 1.5 times the maximum took it above in 27%, by up
to 0.115.

