using Vasilek

include(joinpath(@__DIR__, "scheme_cases.jl"))

"""
Order-of-accuracy suite.

The other advection tests pin a single resolution against a magic `atol`. That
catches regressions but says nothing about whether a scheme has the order it
claims, so it cannot tell a correct scheme from a nearly correct one.

A smooth periodic profile is advected exactly one domain length at a *fixed*
Courant number, so Δt ∝ Δx and the exact solution is the initial condition.
`c = 0.4` keeps the step count integral at every resolution used.

The orders below are **global**, and sit one below the per-step order: the
local error is O(Δxᵖ⁺¹) and the number of steps grows as 1/Δx. Every value was
measured before it was asserted.
"""
const CONV_C = 0.4
const CONV_NS = [32, 64, 128, 256]

conv_norms(e, Δx) = (Δx*sum(abs, e), sqrt(Δx*sum(abs2, e)), maximum(abs, e))

function conv_error(scheme, N)
    Δx = 1/N
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]
    final = march!(f₀, scheme, CONV_C, round(Int, N/CONV_C))
    return conv_norms(final .- f₀, Δx)
end

"Least-squares slope of log(error) against log(Δx) in the given norm."
function conv_order(scheme, which)
    errs = [conv_error(scheme, N)[which] for N in CONV_NS]
    X = log.([1/N for N in CONV_NS])
    y = log.(errs)
    n = length(X)
    return (n*sum(X.*y) - sum(X)*sum(y)) / (n*sum(X.^2) - sum(X)^2)
end

function check_order(name, scheme, expected; which = 2, atol = 0.15)
    p = conv_order(scheme, which)
    nrm = which == 1 ? "L1" : which == 2 ? "L2" : "Linf"
    println("  ", rpad(name, 26), " order = ", lpad(round(p; digits = 2), 5),
            "  (expected ", expected, ", ", nrm, ")")
    @test isapprox(p, expected; atol = atol)
end

@testset "Order of accuracy" begin
    check_order("Upwind", Upwind(), 1)
    check_order("LaxWendroff", LaxWendroff(), 2)
    check_order("Godunov constant", Godunov(PiecewiseConstant()), 1)

    # The limiter clips smooth extrema, so L∞ degrades further (measured 0.69).
    # L1 is the honest norm to hold it to.
    check_order("Godunov linear VanLeer", Godunov(PiecewiseLinear(), VanLeer()), 1; which = 1)

    check_order("SemiLagrangian linear", SemiLagrangian(LinearSpline()), 1)
    check_order("SemiLagrangian quadratic", SemiLagrangian(QuadraticSpline()), 2)
    check_order("SemiLagrangian cubic", SemiLagrangian(CubicSpline()), 3)
    check_order("PFC", PFC(fmin = 0.0, fmax = 2.0), 3)
end

@testset "Scheme equivalences" begin
    # Algebraic identities, not coincidences, and they hold to one ULP.
    N = 64; Δx = 1/N; c = 0.4
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]
    onestep(scheme) = march!(f₀, scheme, c, 1)

    upwind = onestep(Upwind())

    # Godunov with piecewise-constant reconstruction *is* upwind.
    @test onestep(Godunov(PiecewiseConstant())) ≈ upwind rtol=1e-15

    # For 0 < c < 1, linear interpolation at x - cΔx gives
    # f_i(1-c) + f_{i-1}c = f_i - c(f_i - f_{i-1}), which is the upwind formula.
    @test onestep(SemiLagrangian(LinearSpline())) ≈ upwind rtol=1e-15
end

@testset "PiecewiseLinear without a limiter is unconditionally unstable" begin
    # With NoLimiter the limiter returns 1.0 identically, so the flux collapses
    # to |c|*(f[i-1] + f[i])/2 -- a centred flux with forward Euler in time.
    # Von Neumann gives |1 - i·c·sin(kΔx)| = √(1 + c²sin²(kΔx)) > 1 for every
    # mode: unconditionally unstable, worst at kΔx = π/2 where the per-step
    # growth is √(1 + 0.16) ≈ 1.077.
    #
    # A single-step test cannot reveal this. At N = 512 the error reaches
    # 6.66e+24 after a full domain traversal.
    #
    # Documented rather than fixed: stabilising it needs a different time
    # integrator, which is a feature, not a repair. This is why the type carries
    # a warning in its docstring.
    N = 128; Δx = 1/N; c = 0.4
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]

    unlimited = march!(f₀, Godunov(PiecewiseLinear()), c, round(Int, 4N/c))
    growth = maximum(abs, unlimited)/maximum(abs, f₀)
    println("  Godunov PiecewiseLinear, no limiter: amplitude growth ", growth)
    @test growth > 10

    limited = march!(f₀, Godunov(PiecewiseLinear(), VanLeer()), c, round(Int, 4N/c))
    @test maximum(abs, limited) ≤ maximum(abs, f₀) + 1e-12
end
