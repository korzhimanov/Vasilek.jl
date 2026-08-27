using Vasilek

include(joinpath(@__DIR__, "..", "scheme_cases.jl"))

"""
    check_step(name, scheme, Δx, Δt, v, f₀, f₁, tol)

One step of `scheme` on `f₀`, compared with the expected `f₁`.

Single-resolution checks against hand-chosen tolerances. They pin behaviour;
the order-of-accuracy suite is what establishes that a scheme is what it claims
to be. Both signs of `v` are exercised throughout, because `c > 0` and `c < 0`
are separate branches in every scheme.
"""
function check_step(name, scheme, Δx, Δt, v, f₀, f₁, tol)
    dev = Δx*norm(march!(f₀, scheme, v*Δt/Δx, 1) - f₁)
    println("  ", rpad(name, 26), dev)
    @test dev ≈ 0 atol=tol
end

"Schemes and their tolerances for one dataset; `fmax` bounds the PFC limiter."
function scheme_table(t, fmax)
    return [("Upwind",                   Upwind(),                                t.upwind),
            ("LaxWendroff",              LaxWendroff(),                           t.lw),
            ("Godunov constant",         Godunov(PiecewiseConstant()),            t.gc),
            ("Godunov linear",           Godunov(PiecewiseLinear()),              t.gl),
            ("Godunov VanLeer",          Godunov(PiecewiseLinear(), VanLeer()),   t.glv),
            ("SemiLagrangian linear",    SemiLagrangian(LinearSpline()),          t.sl1),
            ("SemiLagrangian quadratic", SemiLagrangian(QuadraticSpline()),       t.sl2),
            ("SemiLagrangian cubic",     SemiLagrangian(CubicSpline()),           t.sl3),
            ("PFC",                      PFC(fmin = 0.0, fmax = fmax),            t.pfc)]
end

@testset "Test 1D advection solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    v = 1.0

    smooth = (
        ("sine",
         [1.0 + 0.01*sin(2π*i*Δx) for i = 0:99],
         [1.0 + 0.01*sin(2π*(i*Δx - v*Δt)) for i = 0:99],
         (upwind = 3e-7, lw = 1e-8, gc = 1e-6, gl = 1e-6, glv = 1e-6,
          sl1 = 3e-7, sl2 = 2e-9, sl3 = 2e-11, pfc = 3e-8)),
        ("gaussian",
         [exp(-((i*Δx - 0.5)/0.15)^2) for i = 0:100],
         [exp(-((i*Δx - v*Δt - 0.5)/0.15)^2) for i = 0:100],
         (upwind = 1e-4, lw = 1e-5, gc = 1e-4, gl = 1.1e-4, glv = 1.2e-4,
          sl1 = 3e-5, sl2 = 4e-7, sl3 = 4e-8, pfc = 5e-6)),
    )

    for (label, f₀, f₁, tols) in smooth
        @testset "$label, forward" begin
            for (name, scheme, tol) in scheme_table(tols, maximum(f₀))
                check_step(name, scheme, Δx, Δt, v, f₀, f₁, tol)
            end
        end
        @testset "$label, backward" begin
            for (name, scheme, tol) in scheme_table(tols, maximum(f₁))
                check_step(name, scheme, Δx, Δt, -v, f₁, f₀, tol)
            end
        end
    end

    # A square pulse at c = 0.8. The monotone schemes reproduce the exact
    # donor-cell answer; the non-monotone ones must not be held to it.
    f₀ = zeros(Float64, 100); f₀[40:50] .= 1.0
    f₁ = zeros(Float64, 100); f₁[41:50] .= 1.0; f₁[40] = 0.2; f₁[51] = 0.8
    f₂ = zeros(Float64, 100); f₂[40:49] .= 1.0; f₂[39] = 0.8; f₂[50] = 0.2

    pulse_tols = (upwind = 1e-18, lw = 1.6e-3, gc = 1e-18, gl = 1e-2, glv = 1e-2,
                  sl1 = 5e-17, sl2 = 2e-3, sl3 = 2e-3, pfc = 1e-18)

    @testset "square pulse, forward" begin
        for (name, scheme, tol) in scheme_table(pulse_tols, 1.0)
            check_step(name, scheme, Δx, Δt, v, f₀, f₁, tol)
        end
    end
    @testset "square pulse, backward" begin
        for (name, scheme, tol) in scheme_table(pulse_tols, 1.0)
            check_step(name, scheme, Δx, Δt, -v, f₀, f₂, tol)
        end
    end
end

@testset "a scheme value is independent of the data it is applied to" begin
    # Under the old closure API the solver captured the arrays it was built
    # with, and two schemes took their loop bounds from the captured array
    # rather than the one passed -- so a differently sized pair read and wrote
    # out of range. Nothing exercised it.
    #
    # A scheme value now holds no data arrays at all, which makes that class of
    # bug unrepresentable. Asserted directly: one value used on two sizes must
    # agree with fresh values used on each.
    big = [1.0 + 0.5*sin(2π*(i-1)/100) for i = 1:100]
    small = [1.0 + 0.5*sin(2π*(i-1)/50) for i = 1:50]

    for scheme in (LaxWendroff(), Upwind(), PFC(fmin = 0.0, fmax = 2.0),
                   SemiLagrangian(CubicSpline()))
        @test march!(big, scheme, 0.4, 2) == march!(big, deepcopy(scheme), 0.4, 2)
        @test march!(small, scheme, 0.4, 2) == march!(small, deepcopy(scheme), 0.4, 2)
        @test length(march!(small, scheme, 0.4, 2)) == 50
    end
end
