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
    # grid scale, so it surfaces late: measured max|f| is still 1.506 after 200
    # steps at c = 1.05, but 1.5e66 after 2000.
    late = march!(f₀, Upwind(), 1.05, 2000)
    println("  Upwind at c = 1.05 after 2000 steps: max|f| = ", maximum(abs, late))
    @test maximum(abs, late) > 1e10
end
