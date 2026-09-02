using Vasilek

"""
Von Neumann amplification factors, against closed forms.

Every scheme here is linear, so on a periodic grid a pure Fourier mode is an
eigenvector: one step multiplies `exp(ikx)` by a complex number `g(kΔx, c)` and
returns a pure mode of the same wavenumber. That number is available in closed
form, and comparing against it is the sharpest statement about these kernels
that this repository can make.

**Why this is not covered by what is already here.** `test_convergence` fits a
slope through four resolutions and asserts it to ±0.15, which is a two-decimal
statement about the aggregate. `test_golden` is bit-exact but pins one dataset
at one Courant number for eight steps, so it detects change without saying what
is right. This suite pins **every mode from the fundamental to the grid scale,
independently, against an analytic value** -- and the same assertion carries
two facts at once: that the output matches `g` and, because the comparison runs
elementwise over `j`, that the output is a pure mode at all. A scheme that is
not translation-invariant fails the second even where it passes the first.

Measured agreement, worst over `m ∈ {1,2,4,8,16,24,31}` and `c ∈ {0.4, 0.8}` at
`N = 64`: 8.9e-15 for the three-point schemes, 4.6e-14 for `PFC`.

**The negative branch is not repeated here.** `test_symmetry` establishes
`R∘A(+c) == A(−c)∘R` to machine precision, which carries every result below to
`c < 0` without a second set of formulas that could disagree with the first.
"""

const AMP_N = 64
const AMP_MODES = (1, 2, 4, 8, 16, 24, 31)
const AMP_COURANTS = (0.4, 0.8)

"""
    amplification(scheme, N, m, c; mean = 0.0, amplitude = 1.0)

The per-`j` ratio `A(f)ⱼ / fⱼ` for `f = mean + amplitude·exp(2πi·m·(j−1)/N)`,
obtained by advecting the real and imaginary parts separately -- the schemes
take real data, and they are real-linear, so the two recombine.

`mean` exists for `PFC`, which needs its data inside `[fmin, fmax]` and away
from the bounds if the limiter is not to engage; it is a fixed point of every
scheme here (`test_symmetry` asserts exactly that) and so is subtracted back
out before the ratio is formed.
"""
function amplification(scheme, N, m, c; mean = 0.0, amplitude = 1.0)
    θ = 2π*m/N
    re = [mean + amplitude*cos(θ*(j-1)) for j = 1:N]
    im = [mean + amplitude*sin(θ*(j-1)) for j = 1:N]
    dr, di = similar(re), similar(im)
    advect!(dr, re, scheme, c, workspace(scheme, N))
    advect!(di, im, scheme, c, workspace(scheme, N))
    return complex.(dr .- mean, di .- mean) ./ complex.(re .- mean, im .- mean)
end

# ------------------------------------------------------------- closed forms
#
# All for the `c > 0` branch, with `θ = kΔx`. Derived from the update formulas,
# not quoted: each one is the symbol of the stencil in `src/VlasovSolver/schemes/`.

"`fᵢ − c(fᵢ − fᵢ₋₁)`."
g_upwind(θ, c) = 1 - c*(1 - exp(-1im*θ))

"`fᵢ − c/2(fᵢ₊₁ − fᵢ₋₁) + c²/2(fᵢ₊₁ − 2fᵢ + fᵢ₋₁)`."
g_lax_wendroff(θ, c) = 1 - 1im*c*sin(θ) + c^2*(cos(θ) - 1)

"""
`Godunov(PiecewiseLinear(), NoLimiter())`.

With the limiter identically 1 the interface value collapses to `(fᵢ₋₁ + fᵢ)/2`
-- a centred flux -- and the update to `fᵢ + c/2(fᵢ₋₁ − fᵢ₊₁)`, whose symbol is
purely `1 − ic·sinθ`. `|g| = √(1 + c²sin²θ) > 1` for every mode that is not
constant or exactly at the grid scale, which is the unconditional instability
`PiecewiseLinear`'s docstring warns about, stated analytically.
"""
g_godunov_linear(θ, c) = 1 - 1im*c*sin(θ)

"""
`PFC` in the regime where its limiter does not engage.

`_ϵ⁺` and `_ϵ⁻` reduce to the plain differences `f₊ − f₀` and `f₀ − f₋` as long
as the data sits well inside `[fmin, fmax]`, leaving the third-order flux

    Φ = c(f₀ + (1−c)/3·[(f₊−f₀)(2−c)/2 + (f₀−f₋)(1+c)/2])

and `dest[i] = src[i] + Φᵢ₋₁ − Φᵢ`. Agreement with this form is what shows the
limiter is inactive on smooth low-amplitude data, which nothing else asserts.
"""
function g_pfc(θ, c)
    Φ = c*(1 + (1 - c)/3*((exp(1im*θ) - 1)/2*(2 - c) +
                          (1 - exp(-1im*θ))/2*(1 + c)))
    return 1 + (exp(-1im*θ) - 1)*Φ
end

"Schemes with a closed-form symbol, and the data they need."
amp_cases() = [
    ("Upwind",             Upwind(),                              g_upwind,          0.0, 1.0),
    ("Godunov constant",   Godunov(PiecewiseConstant()),          g_upwind,          0.0, 1.0),
    ("SemiLagrangian lin", SemiLagrangian(LinearSpline()),        g_upwind,          0.0, 1.0),
    ("LaxWendroff",        LaxWendroff(),                         g_lax_wendroff,    0.0, 1.0),
    ("Godunov linear",     Godunov(PiecewiseLinear()),            g_godunov_linear,  0.0, 1.0),
    ("PFC",                PFC(fmin = 0.0, fmax = 2.0),           g_pfc,             1.0, 0.01),
]

@testset "Amplification factor against the closed form" begin
    # Godunov(PiecewiseConstant) and SemiLagrangian(LinearSpline) share upwind's
    # symbol. `test_convergence` already asserts both equivalences, but on one
    # dataset after one step; here they hold mode by mode, which is the stronger
    # statement and the one that says *why* they are equal.
    worst = Dict{String,Float64}()
    for (name, scheme, g, mean, amplitude) in amp_cases()
        for c in AMP_COURANTS, m in AMP_MODES
            z = amplification(scheme, AMP_N, m, c; mean = mean, amplitude = amplitude)
            d = maximum(abs, z .- g(2π*m/AMP_N, c))
            worst[name] = max(get(worst, name, 0.0), d)
            @test d < 1e-13
        end
    end
    for (name, _, _, _, _) in amp_cases()
        println("  ", rpad(name, 22), "worst |g_num − g_analytic| = ", worst[name])
    end
end

@testset "Dissipation and dispersion" begin
    # The comparison the benchmark suite cannot make: how much amplitude a
    # scheme loses per step, and how far its phase velocity departs from the
    # exact `exp(-icθ)`. Both fall out of `g` at no cost, and both are exact
    # rather than measured, so they can be asserted rather than reported.
    #
    # Measured at c = 0.4, N = 64 (λ = 64Δx down to λ ≈ 2.7Δx):
    #
    #   scheme            |g| m=1     |g| m=16    phase m=1   phase m=16
    #   Upwind            0.998844    0.721110    0.9998      0.9358
    #   LaxWendroff       0.999998    0.930376    0.9987      0.7073
    #   Godunov linear    1.000768    1.077033    0.9979      0.6056
    #   PFC               0.999998    0.890016    1.0000      0.9755
    #
    # where "phase" is arg(g)/(−cθ), which the exact solution makes 1.
    phase(z, θ, c) = angle(z)/(-c*θ)
    g_of(g, m, c) = g(2π*m/AMP_N, c)

    for c in AMP_COURANTS
        println("  c = ", c, "   |g| and arg(g)/(−cθ)")
        for (name, _, g, _, _) in amp_cases()
            row = join([string(rpad(round(abs(g_of(g, m, c)); digits = 6), 8), "/",
                               lpad(round(phase(g_of(g, m, c), 2π*m/AMP_N, c); digits = 4), 7))
                        for m in (1, 4, 16)], "  ")
            println("    ", rpad(name, 22), row, "   (m = 1, 4, 16)")
        end
    end

    for c in AMP_COURANTS
        for m in AMP_MODES
            θ = 2π*m/AMP_N
            # Stable and dissipative: no mode grows, for the three schemes that
            # are stable at 0 < c < 1.
            @test abs(g_of(g_upwind, m, c))       ≤ 1.0 + 1e-14
            @test abs(g_of(g_lax_wendroff, m, c)) ≤ 1.0 + 1e-14
            @test abs(g_of(g_pfc, m, c))          ≤ 1.0 + 1e-14

            # And the unlimited piecewise-linear flux grows on every one of
            # them. This is the analytic statement of what `test_convergence`
            # demonstrates by marching 4N/c steps and watching the amplitude
            # reach 6.66e+24: worst per-step growth 1.077 at c = 0.4.
            @test abs(g_of(g_godunov_linear, m, c)) > 1.0
        end
    end

    @testset "the well-resolved modes separate the schemes" begin
        # At λ = 64Δx upwind has already lost 1.16e-3 per step where
        # LaxWendroff and PFC lose 2e-6 -- three orders of magnitude, and the
        # whole reason a first-order scheme is unusable for a long run.
        for c in AMP_COURANTS
            @test 1 - abs(g_of(g_upwind, 1, c))       > 5e-4
            @test 1 - abs(g_of(g_lax_wendroff, 1, c)) < 1e-5
            @test 1 - abs(g_of(g_pfc, 1, c))          < 1e-5
        end

        # Phase is where PFC separates from LaxWendroff, and it is the opposite
        # ordering to the dissipation above. At λ = 4Δx and c = 0.4 the relative
        # phase velocity is 0.9755 for PFC against 0.7073 for LaxWendroff: a
        # 2.4% error against 29%. LaxWendroff keeps the amplitude of a marginal
        # mode and puts it in the wrong place; upwind removes it instead.
        θ16 = 2π*16/AMP_N
        @test phase(g_of(g_pfc, 16, 0.4), θ16, 0.4) > 0.97
        @test phase(g_of(g_lax_wendroff, 16, 0.4), θ16, 0.4) < 0.75
    end

    @testset "one full traversal, from the symbol alone" begin
        # `|g|^(N/c)` is the amplitude a mode retains after being carried once
        # around the domain -- the exact quantity `test_convergence` measures by
        # marching. Measured for m = 1 at c = 0.4, N = 64: Upwind 0.831004,
        # LaxWendroff 0.999751, PFC 0.999667, Godunov linear 1.130748.
        #
        # The first of those is the 17% amplitude loss that shows up in the
        # work-precision table as upwind's L² error, arrived at independently.
        c = 0.4
        steps = AMP_N/c
        for (name, _, g, _, _) in amp_cases()
            retained = abs(g_of(g, 1, c))^steps
            println("  ", rpad(name, 22), "amplitude after one traversal = ",
                    round(retained; digits = 6))
        end
        @test abs(g_of(g_upwind, 1, c))^steps < 0.85
        @test abs(g_of(g_lax_wendroff, 1, c))^steps > 0.999
        @test abs(g_of(g_pfc, 1, c))^steps > 0.999
        @test abs(g_of(g_godunov_linear, 1, c))^steps > 1.1
    end
end
