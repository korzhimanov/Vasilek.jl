using Vasilek: Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC

"""
Order-of-accuracy suite.

The rest of the advection tests pin a single resolution against a magic
`atol`. That catches regressions but says nothing about whether a scheme has
the order it claims, so it cannot tell a correct scheme from a nearly correct
one.

A smooth periodic profile is advected exactly one domain length at a *fixed*
Courant number, so Δt ∝ Δx and the exact solution is the initial condition.
`c = 0.4` keeps the step count integral at every resolution used.

The orders below are **global**, and are one lower than the per-step order:
the local error is O(Δxᵖ⁺¹) and the number of steps grows as 1/Δx. Every
value here was measured before it was asserted.
"""
const CONV_C = 0.4
const CONV_NS = [32, 64, 128, 256]

conv_norms(e, Δx) = (Δx*sum(abs, e), sqrt(Δx*sum(abs2, e)), maximum(abs, e))

function conv_error(mk, N)
    Δx = 1/N
    x = [(i-1)*Δx for i = 1:N]
    f₀ = @. 1.0 + 0.5*sin(2π*x)
    src = copy(f₀)
    dst = similar(src)
    step! = mk(src, dst)
    for _ = 1:round(Int, N/CONV_C)
        step!(CONV_C)
        copyto!(src, dst)
    end
    return conv_norms(src .- f₀, Δx)
end

"Least-squares slope of log(error) against log(Δx) in the given norm."
function conv_order(mk, which)
    errs = [conv_error(mk, N)[which] for N in CONV_NS]
    X = log.([1/N for N in CONV_NS])
    y = log.(errs)
    n = length(X)
    return (n*sum(X.*y) - sum(X)*sum(y)) / (n*sum(X.^2) - sum(X)^2)
end

function check_order(name, mk, expected; which = 2, atol = 0.15)
    p = conv_order(mk, which)
    nrm = which == 1 ? "L1" : which == 2 ? "L2" : "Linf"
    println("  ", rpad(name, 26), " order = ", lpad(round(p; digits = 2), 5),
            "  (expected ", expected, ", ", nrm, ")")
    @test isapprox(p, expected; atol = atol)
end

@testset "Order of accuracy" begin
    check_order("Upwind", (a,b) -> Upwind.generate_solver(a,b), 1)
    check_order("LaxWendroff", (a,b) -> LaxWendroff.generate_solver(a,b), 2)
    check_order("Godunov constant", (a,b) -> Godunov.generate_solver(a,b,:Riemann_constant), 1)

    # The limiter clips smooth extrema, so L∞ degrades further (measured 0.69).
    # L1 is the honest norm to hold it to.
    check_order("Godunov linear VanLeer",
                (a,b) -> Godunov.generate_solver(a,b,:Riemann_linear; flux_limiter=:VanLeer),
                1; which = 1)

    check_order("SemiLagrangian linear",
                (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Linear), 1)
    check_order("SemiLagrangian quadratic",
                (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Quadratic), 2)
    check_order("SemiLagrangian cubic",
                (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Cubic), 3)
    check_order("PFC", (a,b) -> PFC.generate_solver(a,b;fₘᵢₙ=0.0,fₘₐₓ=2.0), 3)
end

@testset "Scheme equivalences" begin
    # These are algebraic identities, not coincidences, and they hold to one ULP.
    N = 64; Δx = 1/N; c = 0.4
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]
    onestep(mk) = (src = copy(f₀); dst = similar(src); mk(src, dst)(c); dst)

    upwind = onestep((a,b) -> Upwind.generate_solver(a,b))

    # Godunov with piecewise-constant reconstruction *is* upwind.
    @test onestep((a,b) -> Godunov.generate_solver(a,b,:Riemann_constant)) ≈ upwind rtol=1e-15

    # For 0 < c < 1, linear interpolation at x - cΔx gives
    # f_i(1-c) + f_{i-1}c = f_i - c(f_i - f_{i-1}), which is the upwind formula.
    @test onestep((a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Linear)) ≈ upwind rtol=1e-15
end

@testset "Godunov :Riemann_linear is unconditionally unstable" begin
    # With flux_limiter = :None the limiter function returns 1.0 identically, so
    # the flux collapses to |c|*(f[i-1] + f[i])/2 -- a centred flux with forward
    # Euler in time. Von Neumann gives |1 - i·c·sin(kΔx)| = √(1 + c²sin²(kΔx)) > 1
    # for every mode: unconditionally unstable, worst at kΔx = π/2 where the
    # per-step growth is √(1 + 0.16) ≈ 1.077.
    #
    # The existing single-step test passes because one step cannot reveal this.
    # Round-off seeds the grid-scale mode and it takes over: at N = 512 the error
    # reaches 6.66e+24 after a full domain traversal.
    #
    # This is documented rather than fixed. Stabilising it needs a different time
    # integrator, which is a feature, not a repair.
    N = 128; Δx = 1/N; c = 0.4
    f₀ = [1.0 + 0.5*sin(2π*(i-1)*Δx) for i = 1:N]
    src = copy(f₀); dst = similar(src)
    step! = Godunov.generate_solver(src, dst, :Riemann_linear)
    for _ = 1:round(Int, 4N/c)
        step!(c)
        copyto!(src, dst)
    end
    growth = maximum(abs, src)/maximum(abs, f₀)
    println("  Godunov linear (no limiter) amplitude growth: $growth")
    @test growth > 10

    # The limited version stays bounded on exactly the same problem.
    src = copy(f₀); dst = similar(src)
    step! = Godunov.generate_solver(src, dst, :Riemann_linear; flux_limiter = :VanLeer)
    for _ = 1:round(Int, 4N/c)
        step!(c)
        copyto!(src, dst)
    end
    @test maximum(abs, src) ≤ maximum(abs, f₀) + 1e-12
end
