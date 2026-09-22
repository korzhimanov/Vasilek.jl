using Vasilek

@isdefined(march!) || include(joinpath(@__DIR__, "scheme_cases.jl"))

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
function conv_order(scheme, which; Ns = CONV_NS)
    errs = [conv_error(scheme, N)[which] for N in Ns]
    X = log.([1/N for N in Ns])
    y = log.(errs)
    n = length(X)
    return (n*sum(X.*y) - sum(X)*sum(y)) / (n*sum(X.^2) - sum(X)^2)
end

function check_order(name, scheme, expected; which = 2, atol = 0.15, Ns = CONV_NS)
    p = conv_order(scheme, which; Ns = Ns)
    nrm = which == 1 ? "L1" : which == 2 ? "L2" : "Linf"
    span = Ns == CONV_NS ? "" : string(", N = ", first(Ns), "..", last(Ns))
    println("  ", rpad(name, 26), " order = ", lpad(round(p; digits = 2), 5),
            "  (expected ", expected, ", ", nrm, span, ")")
    @test isapprox(p, expected; atol = atol)
end

@testset "Order of accuracy" begin
    check_order("Upwind", Upwind(), 1)
    check_order("LaxWendroff", LaxWendroff(), 2)
    check_order("Godunov constant", Godunov(PiecewiseConstant()), 1)

    # Without a limiter it is LaxWendroff (see below), so only the limited
    # scheme's order needs measuring. VanLeer flattens the slope at a smooth
    # extremum, which costs L∞ and L² more than L¹: measured 2.09 in L¹, 1.77
    # in L² and 1.42 in L∞. L¹ is the honest norm to hold it to.
    #
    # Before the flux carried its (1 − |c|) factor the update was forward Euler
    # on the limited slope, first order at a fixed Courant number -- 1.00, 1.01
    # and 0.82 -- and was held here to 1 in L¹. The 0.82 had been quoted as 0.69.
    check_order("Godunov linear VanLeer", Godunov(PiecewiseLinear(), VanLeer()), 2; which = 1)

    # Superbee is second order in L¹ too, but it gets there late: its limiter
    # is the most compressive in Sweby's region, which costs most on a coarse
    # grid. The local slopes from N = 32 to 2048 are 1.19, 1.81, 1.91,
    # 1.96, 1.98 and 1.99, so a fit over the usual 32 to 256 reads 1.65, still
    # pre-asymptotic. Held to 2 over 128 to 1024, where it measures 1.95. L²
    # settles near 1.69 and L∞ near 1.3.
    check_order("Godunov linear Superbee", Godunov(PiecewiseLinear(), Superbee()), 2;
                which = 1, Ns = [128, 256, 512, 1024])

    check_order("SemiLagrangian linear", SemiLagrangian(LinearSpline()), 1)
    check_order("SemiLagrangian quadratic", SemiLagrangian(QuadraticSpline()), 2)
    check_order("SemiLagrangian cubic", SemiLagrangian(CubicSpline()), 3)
    check_order("PFC", PFC(fmin = 0.0, fmax = 2.0), 3)
end

@testset "Scheme equivalences" begin
    # Algebraic identities, not coincidences, and they hold to two ULP:
    # measured one for Godunov constant, two for the other two.
    N = 64; Δx = 1/N; c = 0.4
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]
    onestep(scheme) = march!(f₀, scheme, c, 1)

    upwind = onestep(Upwind())

    # Godunov with piecewise-constant reconstruction *is* upwind.
    @test onestep(Godunov(PiecewiseConstant())) ≈ upwind rtol=1e-15

    # For 0 < c < 1, linear interpolation at x - cΔx gives
    # f_i(1-c) + f_{i-1}c = f_i - c(f_i - f_{i-1}), which is the upwind formula.
    @test onestep(SemiLagrangian(LinearSpline())) ≈ upwind rtol=1e-15

    # And with piecewise-linear reconstruction and no limiter it *is*
    # Lax–Wendroff: the slope is then the downwind difference, and its average
    # over the strip crossing each interface in a step is Lax–Wendroff's flux,
    # c(f_{i-1} + (1-c)(f_i - f_{i-1})/2). Measured two ULP, 2.2e-16: the two
    # evaluate the same polynomial in a different order.
    @test onestep(Godunov(PiecewiseLinear())) ≈ onestep(LaxWendroff()) rtol=1e-15
end

@testset "PiecewiseLinear without a limiter is stable: it is Lax–Wendroff" begin
    # Until the flux carried its (1 − |c|) factor, the limiter's 1.0 collapsed
    # it to |c|*(f[i-1] + f[i])/2, a centred flux with forward Euler in time:
    # |g| = √(1 + c²sin²(kΔx)) > 1 for every mode, and 1.4e25 times the initial
    # amplitude by the end of this run. The one-step equivalence above holds
    # for the whole run: measured 0.999946 of the initial amplitude at the end,
    # and 1.4e-14 from LaxWendroff after all 1280 steps.
    N = 128; Δx = 1/N; c = 0.4
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]
    steps = round(Int, 4N/c)

    unlimited = march!(f₀, Godunov(PiecewiseLinear()), c, steps)
    lw = march!(f₀, LaxWendroff(), c, steps)
    growth = maximum(abs, unlimited)/maximum(abs, f₀)
    println("  Godunov PiecewiseLinear, no limiter: amplitude ", growth,
            ", max|Δ| from LaxWendroff ", maximum(abs, unlimited .- lw))
    @test growth ≤ 1
    @test maximum(abs, unlimited .- lw) ≤ 1e-12

    limited = march!(f₀, Godunov(PiecewiseLinear(), VanLeer()), c, steps)
    @test maximum(abs, limited) ≤ maximum(abs, f₀) + 1e-12
end
