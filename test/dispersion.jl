# The kinetic dispersion relation, solved rather than tabulated.
#
# Two of the physics tests in this suite used to compare a measured rate against
# a constant that had to come from somewhere else. `LANDAU_ROOTS` was a table of
# three roots copied into the test file, and the two-stream case was held against
# the *cold* closed form because the warm one has no elementary solution. Both
# were workarounds for the same missing piece: the plasma dispersion function.
#
# With `Z` in hand the Landau roots are computed at any `k` rather than quoted at
# three, and the two-stream measurement can be compared against the beams the
# solver was actually given -- warm ones -- instead of against a cold limit it is
# known to miss. Measured at the `vt = 0.3` of the three growth-rate runs, the
# warm root takes the worst error from the cold form's 3.14% to 1.95%. The gap is
# widest at the `vt = 0.6` run, 7.44% against 2.01%, because the cold form's error
# grows with the beam temperature and the warm root's does not.
#
# Deliberately a separate file rather than part of `verification_harness.jl`:
# nothing here runs a simulation, and `test_dispersion.jl` needs it without
# paying for the harness's Strang loop.
#
# The files in `test/` include it as `@isdefined(landau_root) || include(...)`,
# so that each still runs on its own while a full test run, which puts them all
# in `Main`, evaluates this one once: a repeat would redefine every method here,
# and `Pkg.test` runs with `--warn-overwrite=yes`. The guard is `landau_root`
# and not `Z`, the first name defined here: a session with a `Z` of its own -- a
# charge number, say -- would skip this file on `Z`, and the tests would run
# against that `Z` instead. Measured on `test_dispersion.jl`: 12 failures and 2
# errors for a function `Z(x) = x`, 14 errors for `Z = 1`. Included, the file
# replaces the first and refuses the second on the spot, saying why.

using SpecialFunctions: erfcx

"""
    Z(ζ)

The plasma dispersion function of Fried and Conte, `Z(ζ) = i√π·w(ζ)`.

`w` is the Faddeeva function, which `erfcx(-iζ)` gives directly. That matters
rather more than it looks: `Z` is *defined* by a Cauchy integral along the real
axis, which converges only for `Im ζ > 0`, and every root this file goes looking
for is a damped or growing mode with `Im ζ ≠ 0`. `w` is entire, so the
expression below **is** the analytic continuation rather than an approximation
to it, and one formula covers the whole plane.

`erfcx` rather than `erf` because `w(ζ) = exp(-ζ²)erfc(-iζ)` breaks on the real
axis at `|ζ| ≳ 27`: the first factor underflows toward zero, the second
overflows to `Inf`, and `w` comes back with an infinite imaginary part where the
scaled form stays finite. Below that the two agree to the last bit -- measured
identical at `ζ = 8`, 20 and 26 -- so the choice is about the far tail only. The
resonance at `v = ω/k` sits at `ζ = 2.8` for the cases here, but the fluid limit
of a collisional run pushes it out to ten times that.
"""
Z(ζ) = im*sqrt(π)*erfcx(-im*ζ)

"""
    Zprime(ζ)

`dZ/dζ = -2(1 + ζZ(ζ))`.

An identity, not a numerical derivative -- it follows from `Z' = -2(1 + ζZ)`,
which is what makes the susceptibility below a single call to `Z`.
`test_dispersion.jl` checks it against a finite difference anyway, because it is
the one line here that a typo would leave plausible.
"""
Zprime(ζ) = -2*(1 + ζ*Z(ζ))

"""
    susceptibility(ω, k; density = 1.0, drift = 0.0, vt = 1.0)

Electrostatic susceptibility of one drifting Maxwellian species,

    χ = n/(k²vt²)·(1 + ζZ(ζ)),      ζ = (ω - k·u)/(√2·k·vt)

in this package's units, where the thermal velocity is 1 and `ωₚ = 1` for a
species of unit density. `1 + ζZ(ζ)` is `-Z'(ζ)/2`; it is written out rather
than routed through [`Zprime`](@ref) so that the formula in the docstring is the
formula in the code.

Several species add: `dielectric` sums them.
"""
function susceptibility(ω, k; density = 1.0, drift = 0.0, vt = 1.0)
    ζ = (ω - k*drift)/(sqrt(2)*k*vt)
    return density/(k^2*vt^2)*(1 + ζ*Z(ζ))
end

"""
    dielectric(ω, k, species; Δx = 0.0)

`ε(ω, k) = 1 + Σχ`, with `species` a collection of NamedTuples carrying
`density`, `drift` and `vt`. Its roots are the electrostatic modes.

**With `Δx`, it is the dielectric function of the harness's grid rather than of
the continuum.** `make_poisson` solves for the potential spectrally, which is
exact, and then takes the field from it by a centred difference, which is not:
a mode `exp(ikx)` gets `s = sin(kΔx)/(kΔx)` of the field its charge carries. The
field is what the particles feel, so every susceptibility reaches the relation
scaled by the same factor, `ε = 1 + s·Σχ`. That is the one piece of the grid's
error that can be written into the relation exactly; the advection's
truncation error and the splitting's are left to the run. In the bump-on-tail
runs it is almost the whole difference between the run and the continuum root
(see [`bump_on_tail_root`](@ref)).
"""
dielectric(ω, k, species; Δx = 0.0) =
    1 + centred_field(k, Δx)*sum(susceptibility(ω, k; density = s.density,
                                                drift = s.drift, vt = s.vt)
                                 for s in species)

"`sin(kΔx)/(kΔx)`, the fraction of a mode's field a centred difference returns; 1 at `Δx = 0`."
centred_field(k, Δx) = Δx == 0 ? one(float(k)) : sin(k*Δx)/(k*Δx)

"""
    kinetic_root(D, z₀, z₁; tol = 1e-12, maxiter = 100)

Root of the complex function `D` by the secant method, from two starting points.

Secant rather than Newton because `D` is `Z` composed with arithmetic and its
derivative is another `Z` evaluation away; the secant converges superlinearly
here (measured: 6 to 9 iterations from the guesses this file uses) and cannot be
written down wrong. It is not globally convergent, so every caller below passes
a guess it can defend.
"""
function kinetic_root(D, z₀, z₁; tol = 1e-12, maxiter = 100)
    d₀, d₁ = D(z₀), D(z₁)
    for _ in 1:maxiter
        d₁ == d₀ && error("kinetic_root: the secant stalled at $z₁ (D = $d₁)")
        z₂ = z₁ - d₁*(z₁ - z₀)/(d₁ - d₀)
        abs(z₂ - z₁) < tol && return z₂
        z₀, d₀, z₁, d₁ = z₁, d₁, z₂, D(z₂)
    end
    error("kinetic_root: no convergence in $maxiter iterations, last step " *
          "$(abs(z₁ - z₀)) at $z₁")
end

"""
    landau_root(k)

The least damped root of `1 + χ = 0` for a single unit Maxwellian: `ω - iγ`.

This is what `LANDAU_ROOTS` in `test_verification.jl` used to carry as three
typed-in constants. Measured against that table -- the values every Landau
assertion in the suite is held to -- it agrees to every digit the table quotes:

    k     computed                     table
    0.3   1.159846 - 0.012620im        1.15985, 0.01262
    0.4   1.285057 - 0.066128im        1.28506, 0.06613
    0.5   1.415662 - 0.153359im        1.41566, 0.15336

The guess is Bohm--Gross with a nudge into the lower half plane, which is where
the root is for every `k` this package runs: the damping is what the guess does
not know, and the secant finds it in under ten steps.
"""
landau_root(k) = kinetic_root(ω -> dielectric(ω, k, ((density = 1.0, drift = 0.0, vt = 1.0),)),
                              sqrt(1 + 3k^2) - 0.05im,
                              sqrt(1 + 3k^2) - 0.06im)

"""
    two_stream_warm(a; v₀ = 3.0, vt = 0.3)

Growth rate of the two-stream instability for two **warm** Maxwellian beams of
density 1/2 at `±v₀`, at `a = k·v₀`. Zero when the configuration is stable, so
it drops into the places `γ_cold` occupied.

The unstable root of the symmetric two-beam system is purely imaginary -- the
two beams contribute complex conjugate susceptibilities at `ω = iγ`, so `ε` is
real on the imaginary axis -- which is why this bisects a real function instead
of calling [`kinetic_root`](@ref). Bisection also gives the stability statement
for free: no sign change in `(0, 1]` means no growing root, and the function
returns 0.0 rather than a value nobody should trust.

**This is the comparison `γ_cold` cannot make.** Measured against the rates the
suite fits, with the beams at `vt = 0.3` the runs actually use:

    a     warm      cold      measured   vs warm   vs cold
    0.4   0.30362   0.30819   0.30173    -0.62%    -2.10%
    0.6   0.34909   0.35339   0.34228    -1.95%    -3.14%
    0.8   0.31201   0.31134   0.31228    +0.09%    +0.30%
    1.0   0.09823   0         0.09512    -3.17%      --

Two things fall out of that table which the cold form actively misleads about.
The first is the tolerance: 2.01% against the warm root where the cold one needs
6%. The second is the sign at `a = 0.8`, where the warm rate is *above* the cold
one -- the crossing sits between `a = 0.75` and `0.8`, and beyond `a = 1` the
cold branch is identically zero while warm beams are still unstable, which the
suite already measures independently. "Finite temperature cannot make a beam
grow faster than a cold one" is therefore false near the band edge, and this
function is how that stopped being a matter of opinion.
"""
function two_stream_warm(a; v₀ = 3.0, vt = 0.3)
    k = a/v₀
    beams = ((density = 0.5, drift = v₀, vt = vt), (density = 0.5, drift = -v₀, vt = vt))
    D(γ) = real(dielectric(im*γ, k, beams))
    lo, hi = 1e-6, 1.0
    D(lo)*D(hi) < 0 || return 0.0
    for _ in 1:200
        mid = 0.5*(lo + hi)
        D(lo)*D(mid) ≤ 0 ? (hi = mid) : (lo = mid)
    end
    return 0.5*(lo + hi)
end

"""
    BUMP_ON_TAIL

The bump-on-tail distribution of Arber and Vann [J. Comput. Phys. 180, 339
(2002)], as the two Maxwellians it is made of. The literature writes it as

    F(v) = [0.9·exp(−v²/2) + 0.2·exp(−2(v − 4.5)²)]/√(2π)

which is a bulk of density 0.9 at rest and unit temperature, and a beam drifting
at 4.5 with `vt = 0.5` -- whose density is therefore 0.1, not the 0.2 in front of
it, since a width of 0.5 halves the integral. The two add to 1 exactly.

`test_dispersion.jl` checks this pair against the formula above by a quadrature
that never calls `Z`, and `bump_on_tail` in the harness runs the formula, not
the pair.
"""
const BUMP_ON_TAIL = ((density = 0.9, drift = 0.0, vt = 1.0),
                      (density = 0.1, drift = 4.5, vt = 0.5))

"""
    bump_on_tail_root(k; Δx = 0.0)

The growing root `ω + iγ` of the [`BUMP_ON_TAIL`](@ref) dispersion relation: a
wave travelling with the beam, its phase velocity on the beam's rising flank.

At the `k = 0.3` the runs use,

    ω = 1.001218,  γ = 0.198098,  v_φ = ω/k = 3.337

near the top of the unstable band: the rate peaks at 0.2028 at `k = 0.266`, and
the band closes at `k = 0.4824`, where the phase velocity has come down to the
bottom of the valley between bulk and beam, 3.104 -- the minimum of `F`, as
Penrose's criterion requires of a marginal mode. The box's harmonics, from 0.6
up, are outside it.

**This is the beam-plasma instability, and the beam's temperature corrects it
rather than drives it.** A cold beam on the same warm bulk grows at 0.2333, and
the root cooled continuously from `vt = 0.5` meets that; `γ` here is 1.3 times
`k·v_t` of the beam, far from the weak-beam limit where the growth is Landau
damping run in reverse. The beam's warmth takes 15% off the cold rate, and that
is what a run has to resolve to land on this root rather than the fluid one.

**The guess is the steepest point of the flank**, `ω = k(u_b − v_t,b)` -- where
a Gaussian beam's slope is largest, 4.0 here -- nudged into the upper half plane
by 0.1i. The phase velocity stays on the flank across the band, 3.1 to 3.9, so
that is where the root is for every `k`: from it the secant finds the growing
root at all fourteen `k` from 0.125 to 0.45 in steps of 0.025. A Bohm--Gross guess,
which `landau_root` can afford, misses it for `k ≤ 0.225`, landing on a neutral
real root or on none: the growing wave there oscillates well below `ωₚ` (0.48 at
`k = 0.125`), nowhere near Bohm--Gross.

With `Δx` the root is the grid's (see [`dielectric`](@ref)): at the 64 cells of
the runs, `ω = 1.001052` and `γ = 0.197872`, 0.017% and 0.114% below the
continuum.
"""
function bump_on_tail_root(k; Δx = 0.0)
    beam = BUMP_ON_TAIL[2]
    guess = k*(beam.drift - beam.vt) + 0.1im
    return kinetic_root(ω -> dielectric(ω, k, BUMP_ON_TAIL; Δx = Δx), guess, guess + 0.01)
end
