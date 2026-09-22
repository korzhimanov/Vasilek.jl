using Vasilek

include(joinpath(@__DIR__, "scheme_cases.jl"))

"""
Structural properties of the advection schemes: conservation, monotonicity,
positivity and the stability limit. None of this was covered before.

Every threshold was measured first and is quoted beside it, so a future change
that shifts one is visible rather than silently re-tuned away.
"""
const INV_N = 128
const INV_C = 0.4

inv_smooth() = [1.0 + 0.5*sin(2π*(i-1)/INV_N) for i = 1:INV_N]
inv_pulse() = [0.4 < (i-1)/INV_N < 0.6 ? 1.0 : 0.0 for i = 1:INV_N]

total_variation(f) = sum(abs(f[mod1(i+1, length(f))] - f[i]) for i in eachindex(f))

const INV_SCHEMES = [
    ("Upwind",                Upwind()),
    ("LaxWendroff",           LaxWendroff()),
    ("Godunov constant",      Godunov(PiecewiseConstant())),
    ("Godunov VanLeer",       Godunov(PiecewiseLinear(), VanLeer())),
    ("Godunov Superbee",      Godunov(PiecewiseLinear(), Superbee())),
    ("SemiLagrangian linear", SemiLagrangian(LinearSpline())),
    ("SemiLagrangian cubic",  SemiLagrangian(CubicSpline())),
]

@testset "Mass conservation" begin
    # At constant velocity on a periodic uniform grid every scheme here
    # conserves Σf to round-off after 500 steps. The flux-form schemes
    # (Upwind, Godunov, PFC) do so by construction; LaxWendroff telescopes; and
    # semi-Lagrangian interpolation weights form a partition of unity, so a
    # constant shift preserves the sum too.
    #
    # Measured relative drift: 1.1e-16 (Upwind) through 1.5e-14 (SL cubic).
    #
    # This says nothing about a *varying* velocity, where semi-Lagrangian is not
    # conservative. That case needs its own test and does not have one.
    f₀ = inv_smooth()
    for (name, scheme) in vcat(INV_SCHEMES, [("PFC", PFC(fmin = 0.0, fmax = 2.0))])
        f = march!(f₀, scheme, INV_C, 500)
        drift = abs(sum(f) - sum(f₀))/sum(f₀)
        println("  mass drift ", rpad(name, 24), drift)
        @test drift < 1e-13
    end
end

@testset "Total variation" begin
    # Godunov's theorem: no linear scheme of order above one can be monotone.
    # LaxWendroff and cubic semi-Lagrangian are exactly that, and both overshoot
    # on a square pulse -- measured TV growth 1.87 and 1.28, undershooting to
    # -0.227 and -0.062.
    f₀ = inv_pulse()
    tv₀ = total_variation(f₀)

    for (name, scheme) in INV_SCHEMES
        f = march!(f₀, scheme, INV_C, 200)
        ratio = total_variation(f)/tv₀
        println("  TV ratio ", rpad(name, 24), round(ratio; digits = 4),
                "   min = ", round(minimum(f); digits = 4))
        if name in ("LaxWendroff", "SemiLagrangian cubic")
            @test ratio > 1               # not TVD, and must stay documented as such
            @test minimum(f) < -1e-3      # it genuinely undershoots
        else
            @test ratio ≤ 1 + 1e-12       # TVD
        end
    end

    f = march!(f₀, PFC(fmin = 0.0, fmax = 1.0), INV_C, 200)
    @test total_variation(f)/tv₀ ≤ 1 + 1e-12
end

@testset "Godunov VanLeer and Superbee are TVD up to the Courant limit" begin
    # The testset above runs them at c = 0.4 only. Sweby's flux-limited
    # Lax–Wendroff with a limiter inside his region -- φ(r) ≤ 2 and φ(r)/r ≤ 2,
    # which VanLeer keeps and Superbee reaches -- is TVD for every |c| ≤ 1, and
    # so obeys a discrete maximum principle: each new value is a convex
    # combination of two old ones. Measured over 200 steps of the pulse, both
    # directions, at every c below: TV ratio 0.999992 to 1 for VanLeer and
    # 0.99999999 to 1 for Superbee, with f inside [0, 1] for both; and a sine
    # that touches zero stays above it at every one of 1000 steps.
    #
    # Godunov VanLeer's flux used to lack the (1 − |c|) factor, and the update
    # was then TVD only to |c| ≤ 1/2: over the same 200 steps the TV ratio was
    # 2.05 at 0.6, 39 at 0.7, 4.8e5 at 0.8 and 1.0e7 at 0.9, and the sine went
    # negative in one step from c = 0.58, to -9.0e-6.
    f₀ = inv_pulse()
    tv₀ = total_variation(f₀)
    touching = [0.5*(1 + sin(2π*(i-1)/INV_N)) for i = 1:INV_N]    # min is 0
    for (name, limiter) in (("VanLeer", VanLeer()), ("Superbee", Superbee()))
        scheme = Godunov(PiecewiseLinear(), limiter)
        worst = (ratio = 0.0, min = Inf, max = -Inf, lowest = Inf)
        for c in (0.5, 0.58, 0.6, 0.7, 0.8, 0.9, 1.0), s in (1, -1)
            f = march!(f₀, scheme, s*c, 200)
            ratio = total_variation(f)/tv₀

            g = copy(touching)
            h = similar(g)
            lowest = minimum(g)
            for _ = 1:1000
                advect!(h, g, scheme, s*c)
                g, h = h, g
                lowest = min(lowest, minimum(g))
            end
            worst = (ratio = max(worst.ratio, ratio), min = min(worst.min, minimum(f)),
                     max = max(worst.max, maximum(f)), lowest = min(worst.lowest, lowest))
            @test ratio ≤ 1 + 1e-12
            @test minimum(f) ≥ 0.0
            @test maximum(f) ≤ 1.0
            @test lowest ≥ 0.0
        end
        println("  Godunov ", rpad(name, 9), "c = ±0.5 to ±1: worst TV ratio ", worst.ratio,
                ", f in [", worst.min, ", ", worst.max, "], sine's lowest ", worst.lowest)
    end
end

@testset "Superbee is anti-diffusive on smooth data" begin
    # Compression has a price the TVD property does not show. Superbee steepens
    # a smooth slope, and the L² norm of the perturbation *grows*, where every
    # other scheme here loses some. After one traversal at c = 0.4 and N = 128,
    # it is +1.3e-3 on the sine and +1.2e-2 on `test_comparison`'s gaussian,
    # against -1.1e-4 and -6.9e-3 for VanLeer. The growth shrinks with the grid
    # (+8.9e-5 and +1.4e-3 at N = 512) but has the same sign at N = 64, 128 and
    # 512.
    #
    # In a Vlasov–Poisson run the same growth shows in f's L², and it biases a
    # damping rate low. On the Landau refinement ladder of `test_verification`
    # (k = 0.5), Superbee's error in γ goes -2.41%, -0.39%, -0.68%, +0.13% with
    # L² rising at every level, where VanLeer's falls monotonically, +7.76% to
    # +0.43%.
    l2(f) = sqrt(sum(abs2, f)/length(f))
    for (label, f₀) in (("sine", inv_smooth()),
                        ("gaussian", [1.0 + exp(-((((i-1)/INV_N) - 0.5)/0.08)^2) for i = 1:INV_N]))
        mean = sum(f₀)/INV_N
        change(scheme) = l2(march!(f₀, scheme, INV_C, round(Int, INV_N/INV_C)) .- mean)/
                         l2(f₀ .- mean) - 1
        sb = change(Godunov(PiecewiseLinear(), Superbee()))
        vl = change(Godunov(PiecewiseLinear(), VanLeer()))
        println("  L² of the perturbation after a traversal, ", rpad(label, 9),
                "Superbee ", sb, ", VanLeer ", vl)
        @test sb > 0
        @test vl < 0
    end
end

@testset "PFC positivity and maximum principle" begin
    # Preserving positivity is the entire reason this scheme exists, and there
    # was no test for it. 1000 steps on a square pulse: measured min = 1.4e-31,
    # max = 0.9999994.
    f = march!(inv_pulse(), PFC(fmin = 0.0, fmax = 1.0), INV_C, 1000)
    println("  PFC after 1000 steps: min = ", minimum(f), "  max = ", maximum(f))
    @test minimum(f) ≥ 0.0
    @test maximum(f) ≤ 1.0 + eps()
end

@testset "Courant limit" begin
    f₀ = inv_smooth()

    # c = 1 is a pure one-cell translation: f[i] - 1*(f[i] - f[i-1]) = f[i-1].
    # Exactly, in floating point.
    @test march!(f₀, Upwind(), 1.0, 1) == circshift(f₀, 1)
    # and after a full lap it returns to the initial data bit-for-bit
    @test march!(f₀, Upwind(), 1.0, INV_N) == f₀

    # Beyond it, upwind is unstable. The growth is seeded by round-off at the
    # grid scale -- |1 − 2c| = 1.1 per step for the grid-scale mode at c = 1.05
    # -- so it surfaces late: measured max|f| was still 1.506 after 200 steps,
    # but 1.5e66 after 2000. That lateness is why `advect!` now refuses the step
    # rather than leaving it to be noticed; `test_contracts.jl` covers the
    # refusal for every scheme.
    @test_throws DomainError march!(f₀, Upwind(), 1.05, 1)
end
