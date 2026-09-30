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
# warm root takes the worst error from the cold form's 3.15% to 1.95%. The gap is
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
using LinearAlgebra: det

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
    0.4   0.30362   0.30819   0.30172    -0.63%    -2.10%
    0.6   0.34909   0.35339   0.34227    -1.95%    -3.15%
    0.8   0.31201   0.31134   0.31228    +0.09%    +0.30%
    1.0   0.09823   0         0.09513    -3.16%      --

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

# ------------------------------------------------------------------ collisions
#
# The same relation with the `BGK` operator in the loop. Everything above is
# collisionless; the one collision operator the package exports had been
# checked only on a single velocity line, never in a run with a field.

"""
    resolvent_moments(ζ)

`Jₘ(ζ) = ∫vᵐF₀(v)/(v − ζ)dv` for `m = 0, …, 4`, with `F₀` the unit Maxwellian,
on the Landau contour: the plain integral for `Im ζ > 0`, its continuation
below.

`J₀ = Z(ζ/√2)/√2`, and the rest follow from `v/(v − ζ) = 1 + ζ/(v − ζ)`, which
gives `Jₘ₊₁ = μₘ + ζJₘ` with `μₘ` the Maxwellian's own moments 1, 0, 1, 0. So
`J₁ = 1 + ζJ₀` is the `1 + ζZ` of the collisionless susceptibility, and the
higher ones are what a collision operator that restores more than the density
needs.
"""
function resolvent_moments(ζ)
    J₀ = Z(ζ/sqrt(2))/sqrt(2)
    J₁ = 1 + ζ*J₀
    J₂ = ζ*J₁
    J₃ = 1 + ζ*J₂
    J₄ = ζ*J₃
    return (J₀, J₁, J₂, J₃, J₄)
end

"""
    collisional_determinant(ω, k, ν; s = 1.0, conserve = (:n, :u, :T))

The linear response of a unit Maxwellian to Vlasov--Poisson with the `BGK`
operator at rate `ν = 1/τ`, as a determinant whose zeros are the modes.

`BGK` relaxes each velocity line towards the Maxwellian with the line's own
density, mean velocity and temperature. Linearised about `F₀`, that Maxwellian is

    M₁ = F₀·[n₁ + u₁v + T₁(v² − 1)/2],
    n₁ = ∫f₁dv,   u₁ = ∫vf₁dv,   T₁ = ∫(v² − 1)f₁dv

and a mode `exp(i(kx − ωt))` of the kinetic equation, with the field the
harness's -- `ik·e₁ = n₁`, the kick `+e∂ᵥ` -- is

    f₁ = F₀·[−(is/k)·n₁·v + ν·M₁/F₀] / (ik(v − ζ)),     ζ = (ω + iν)/k

`s` scales the field, as `centred_field` does for [`dielectric`](@ref). Taking the
three moments of `f₁` closes a 3×3 linear system on `(n₁, u₁, T₁)` through the
[`resolvent_moments`](@ref), and a mode is where its determinant vanishes. At
`ν = 0` the last two rows are the identity and the determinant **is**
`1 + (s/k²)(1 + ζZ)`, the collisionless dielectric function.

`conserve` names what the relaxation restores, for the two operators `BGK` is
not: `(:n,)` is the Krook model, relaxing to the background temperature at the
local density, and `(:n, :u)` restores momentum but not energy. They are here
because they are what a slip in the operator would turn it into, and the modes
tell the three apart; see [`collisional_root`](@ref).
"""
function collisional_determinant(ω, k, ν; s = 1.0, conserve = (:n, :u, :T))
    issubset(conserve, (:n, :u, :T)) ||
        throw(ArgumentError("conserve takes :n, :u and :T, got $conserve"))
    ζ = (ω + im*ν)/k
    J₀, J₁, J₂, J₃, J₄ = resolvent_moments(ζ)
    cn, cu, cT = (:n in conserve), (:u in conserve), (:T in conserve)
    # f₁ = F₀·(a₀ + a₁v + a₂v²)/(ik(v − ζ)), each aᵢ linear in (n₁, u₁, T₁):
    # the tuples are their coefficients on the three.
    a₀ = (ν*cn, 0.0, -ν*cT/2)
    a₁ = (-im*s/k, ν*cu, 0.0)
    a₂ = (0.0, 0.0, ν*cT/2)
    A = [(j == i) - (a₀[j]*m₀ + a₁[j]*m₁ + a₂[j]*m₂)/(im*k)
         for (i, (m₀, m₁, m₂)) in enumerate(((J₀, J₁, J₂),
                                             (J₁, J₂, J₃),
                                             (J₂ - J₀, J₃ - J₁, J₄ - J₂))),
             j in 1:3]
    return det(A)
end

"""
    collisional_susceptibility(ω, k, ν; conserve = (:n, :u, :T))

`χ` such that the dielectric function with collisions is `ε = 1 + s·χ`, the form
[`dielectric`](@ref) has: the field enters [`collisional_determinant`](@ref)
through one coefficient, so the determinant is linear in `s`, `D(s) = D₀ + s·D₁`,
and `ε = D(s)/D₀`. `D₀` is the response with the field switched off.

At `ν = 0` this is [`susceptibility`](@ref) exactly.
"""
function collisional_susceptibility(ω, k, ν; conserve = (:n, :u, :T))
    D₀ = collisional_determinant(ω, k, ν; s = 0.0, conserve)
    return (collisional_determinant(ω, k, ν; s = 1.0, conserve) - D₀)/D₀
end

"""
    collisional_root(k, ν; Δx = 0.0, conserve = (:n, :u, :T))

The Langmuir root `ω − iγ` at collision rate `ν`, followed from the Landau root
at `ν = 0` in small steps of `ν`. With `Δx` it is the grid's, as for
[`dielectric`](@ref).

At `k = 0.5`, for the operator `BGK` is and for the two it is not:

    ν      (:n, :u, :T) -- BGK    (:n, :u)               (:n,) -- Krook
    0      1.415662 - 0.153359i   (the Landau root, for all three)
    0.1    1.401930 - 0.132698i   1.380107 - 0.162038i   1.395693 - 0.221399i
    0.3    1.382899 - 0.106214i   1.321691 - 0.166026i   1.348187 - 0.355524i
    1      1.351750 - 0.065176i   1.212495 - 0.137218i   1.078340 - 0.809220i

**Collisions that restore all three moments weaken the damping; restoring fewer
strengthens it.** `dγ/dν` at `ν = 0` is −0.243, +0.114 and +0.684. Landau damping
at `k = 0.5` is resonant particles at `v = 2.83` phase-mixing, and scattering
them spoils the resonance; what replaces it depends on what the collisions keep.
Kept momentum and energy leave a fluid, whose wave damps only by conducting
heat, at a rate that falls as `1/ν`, and the root goes to the adiabatic
`√(1 + 3k²)` (see [`heat_mode_root`](@ref)). Without energy it goes to the
isothermal `√(1 + k²)`, and without momentum collisions are friction on the flow
itself: the Krook pair slows, damps harder, and merges on the imaginary axis
between `ν = 1.8` and 1.9, past which following it means nothing.

So the sign of the change is enough to say which operator a run has, and by
`ν = 1` the three frequencies are 10% and 20% apart.

The zeros are taken of [`collisional_determinant`](@ref) rather than of `ε`,
which has poles where the field-free response `D₀` vanishes.

**The secant stops at steps of 1e-9 rather than [`kinetic_root`](@ref)'s 1e-12**,
because the determinant is not that precise at large `ν`. `J₃` and `J₄` are built
as `1 + ζJ₂`, which cancels to `|ζ|²` once `|ζ| = |ω + iν|/k` is large: measured
against quadrature, they carry 1.2e-10 of their size at `ν = 20` and `k = 0.5`,
and 2.8e-10 at `ν = 40`, where `J₀` to `J₂` carry 3e-13. The 1e-12 steps sit
inside that noise, and at `ν = 20` the secant ran out of iterations taking them.
Superlinear convergence leaves the returned root far better than the last step.
"""
function collisional_root(k, ν; Δx = 0.0, conserve = (:n, :u, :T))
    s = centred_field(k, Δx)
    electrons = ((density = 1.0, drift = 0.0, vt = 1.0),)
    z = kinetic_root(ω -> dielectric(ω, k, electrons; Δx = Δx),
                     sqrt(1 + 3k^2) - 0.05im, sqrt(1 + 3k^2) - 0.06im)
    ν == 0 && return z
    for νⱼ in range(0.0, ν; length = 41 + ceil(Int, 2ν))[2:end]
        z = kinetic_root(ω -> collisional_determinant(ω, k, νⱼ; s, conserve), z, z + 1e-4;
                         tol = 1e-9)
    end
    return z
end

"""
    heat_mode_root(k, ν; Δx = 0.0)

The rate `g` of the purely damped root `ω = −ig` of [`collisional_determinant`](@ref)
with all three moments restored: heat conduction, the one mode here that does
not oscillate.

**It exists because `BGK` conserves energy.** A fluid with an energy equation has
three modes at each `k` -- the two Langmuir waves and a non-propagating
temperature perturbation that decays by conduction -- and one without has two.
Scanned along the axis at `k = 0.5`, `(:n, :u)` has no root at all below `g = 5`
at `ν = 1` or below `g = 14` at `ν = 10`, where BGK's is at 0.434 and 0.0534.

The Chapman--Enskog closure of 1D BGK has no viscosity -- in one dimension the
pressure is `nT` by the definition of `T` -- and conducts heat as
`Q = −(3nT/ν)∂ₓT`, which puts the three modes at the roots of

    ω³ + iχω² − (1 + 3k²)ω − iχ(1 + k²) = 0,      χ = 3k²/ν

so `g → χ(1 + k²)/(1 + 3k²)`, 0.5357/ν at `k = 0.5`, and the Langmuir root
`→ √(1 + 3k²) − iχk²/(1 + 3k²)`: the adiabatic Bohm--Gross frequency, since one
degree of freedom makes the adiabatic index 3, and a damping that falls as
`1/ν`. `test_dispersion.jl` holds the kinetic roots to that limit.

The determinant is real on the imaginary axis -- `F₀` is even, so `D(−ω*) = D(ω)*`
-- which is why this bisects a real function, as [`two_stream_warm`](@ref) does,
from a scan for the first sign change. **The scan stops at `g = ν + 2k`**, where
`Im ζ = −2`: further down, the continuation of `Z` grows as `exp(|ζ|²/2)` and
`1 + ζZ` cancels away every digit, and a scan to `g = 30` at `ν = 1` found
dozens of spurious sign changes past `g ≈ 5`, in BGK's determinant and in that
of `(:n, :u)`. The heat mode is above that line whenever it is a fluid mode at
all.
"""
function heat_mode_root(k, ν; Δx = 0.0)
    s = centred_field(k, Δx)
    D(g) = real(collisional_determinant(-im*g, k, ν; s))
    grid = exp.(range(log(1e-4), log(ν + 2k); length = 400))
    i = findfirst(j -> D(grid[j])*D(grid[j+1]) < 0, 1:length(grid)-1)
    i === nothing && error("heat_mode_root: no sign change of D on the imaginary " *
                           "axis above Im ζ = -2 at ν = $ν")
    lo, hi = grid[i], grid[i+1]
    for _ in 1:200
        mid = 0.5*(lo + hi)
        D(lo)*D(mid) ≤ 0 ? (hi = mid) : (lo = mid)
    end
    return 0.5*(lo + hi)
end
