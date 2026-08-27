using Vasilek: Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC

"""
Structural properties of the advection schemes: conservation, monotonicity,
positivity and the stability limit. None of this was covered before.

Every threshold here was measured first and is quoted in the comment beside
it, so a future change that shifts one is visible rather than silently
re-tuned away.
"""
const INV_N = 128
const INV_C = 0.4

inv_x() = [(i-1)/INV_N for i = 1:INV_N]
inv_smooth() = (x = inv_x(); @. 1.0 + 0.5*sin(2π*x))
inv_pulse() = [0.4 < xi < 0.6 ? 1.0 : 0.0 for xi in inv_x()]

total_variation(f) = sum(abs(f[mod1(i+1, length(f))] - f[i]) for i in eachindex(f))

function inv_march(mk, f₀, nsteps; c = INV_C)
    src = copy(f₀)
    dst = similar(src)
    step! = mk(src, dst)
    for _ = 1:nsteps
        step!(c)
        copyto!(src, dst)
    end
    return src
end

const INV_SCHEMES = [
    ("Upwind",           (a,b) -> Upwind.generate_solver(a,b)),
    ("LaxWendroff",      (a,b) -> LaxWendroff.generate_solver(a,b)),
    ("Godunov constant", (a,b) -> Godunov.generate_solver(a,b,:Riemann_constant)),
    ("Godunov VanLeer",  (a,b) -> Godunov.generate_solver(a,b,:Riemann_linear; flux_limiter=:VanLeer)),
    ("SemiLagrangian linear", (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Linear)),
    ("SemiLagrangian cubic",  (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Cubic)),
]

@testset "Mass conservation" begin
    # At constant velocity on a periodic uniform grid every scheme here
    # conserves Σf to round-off after 500 steps. The flux-form schemes
    # (Upwind, Godunov, PFC) do so by construction; LaxWendroff telescopes;
    # and semi-Lagrangian interpolation weights form a partition of unity, so
    # a constant shift preserves the sum too.
    #
    # Measured relative drift: 6.7e-16 (Godunov constant), 8.9e-16 (VanLeer),
    # 1.2e-15 (PFC), 4.9e-15 (SL linear), 1.5e-14 (SL cubic).
    #
    # This says nothing about a *varying* velocity, where semi-Lagrangian is
    # not conservative. That case needs its own test and does not have one.
    f₀ = inv_smooth()
    for (name, mk) in INV_SCHEMES
        f = inv_march(mk, f₀, 500)
        drift = abs(sum(f) - sum(f₀))/sum(f₀)
        println("  mass drift ", rpad(name, 24), drift)
        @test drift < 1e-13
    end
    f = inv_march((a,b) -> PFC.generate_solver(a,b;fₘᵢₙ=0.0,fₘₐₓ=2.0), f₀, 500)
    @test abs(sum(f) - sum(f₀))/sum(f₀) < 1e-13
end

@testset "Total variation" begin
    # Godunov's theorem: no linear scheme of order above one can be monotone.
    # LaxWendroff and cubic semi-Lagrangian are exactly that, and both
    # overshoot on a square pulse — measured TV growth 1.87 and 1.28, with
    # undershoot to -0.227 and -0.062 respectively.
    f₀ = inv_pulse()
    tv₀ = total_variation(f₀)

    for (name, mk) in INV_SCHEMES
        f = inv_march(mk, f₀, 200)
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

    f = inv_march((a,b) -> PFC.generate_solver(a,b;fₘᵢₙ=0.0,fₘₐₓ=1.0), f₀, 200)
    @test total_variation(f)/tv₀ ≤ 1 + 1e-12
end

@testset "PFC positivity and maximum principle" begin
    # Preserving positivity is the entire reason this scheme exists, and there
    # was no test for it. 1000 steps on a square pulse: measured
    # min = 1.4e-31, max = 0.9999994.
    f₀ = inv_pulse()
    f = inv_march((a,b) -> PFC.generate_solver(a,b;fₘᵢₙ=0.0,fₘₐₓ=1.0), f₀, 1000)
    println("  PFC after 1000 steps: min = ", minimum(f), "  max = ", maximum(f))
    @test minimum(f) ≥ 0.0
    @test maximum(f) ≤ 1.0 + eps()
end

@testset "Courant limit" begin
    f₀ = inv_smooth()

    # c = 1 is a pure one-cell translation: f[i] - 1*(f[i] - f[i-1]) = f[i-1].
    # Exactly, in floating point.
    @test inv_march((a,b) -> Upwind.generate_solver(a,b), f₀, 1; c = 1.0) == circshift(f₀, 1)
    # and after a full lap it returns to the initial data bit-for-bit
    @test inv_march((a,b) -> Upwind.generate_solver(a,b), f₀, INV_N; c = 1.0) == f₀

    # Beyond it, upwind is unstable. The growth is seeded by round-off at the
    # grid scale, so it takes a while to surface: measured max|f| is still
    # 1.506 after 200 steps at c = 1.05, but 1.5e66 after 2000.
    late = inv_march((a,b) -> Upwind.generate_solver(a,b), f₀, 2000; c = 1.05)
    println("  Upwind at c = 1.05 after 2000 steps: max|f| = ", maximum(abs, late))
    @test maximum(abs, late) > 1e10
end
