using Vasilek
using Vasilek.Collisions: Landau1P, collide!, workspace as collision_workspace

include(joinpath(@__DIR__, "scheme_cases.jl"))

"""
What the suite could not answer until now: **how the schemes compare with each
other**, rather than how each one behaves on its own.

`benchmark/` times every scheme in isolation and never relates cost to accuracy.
`test_convergence` establishes that each scheme has the order it claims. Neither
says which scheme a caller should reach for, and the answer is not a single
name -- it inverts between problem classes, which is the thing worth pinning.

**Only deterministic quantities are asserted here.** Errors are reproducible
bit-for-bit on a given machine and across the ones this package is tested on;
wall-clock times are not, at any tolerance worth having, as
`benchmark/runbenchmarks.jl` argues from measurement. So the accuracy *ranking*
is a gate and the cost is not. The one timing statement that is gated is the
complexity class at the bottom, which separates O(N) from O(N²) by a factor of
two per doubling and cannot be crossed by runner noise.

The companion report, which pairs these errors with measured cost, is
`benchmark/workprecision.jl`.
"""

const CMP_C = 0.4
const CMP_N = 512

cmp_sine(N)   = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
cmp_gauss(N)  = [1.0 + exp(-((((i-1)/N) - 0.5)/0.08)^2) for i = 1:N]
cmp_square(N) = [0.4 < (i-1)/N < 0.6 ? 1.5 : 1.0 for i = 1:N]

"""
    traversal_error(scheme, profile, N = CMP_N)

L² error after carrying `profile` exactly once around the periodic domain at a
fixed Courant number, so the exact answer is the initial condition.

The same construction `test_convergence` uses for its order fits; here it is
read across schemes at one resolution rather than across resolutions for one
scheme.
"""
function traversal_error(scheme, profile, N = CMP_N)
    Δx = 1/N
    f₀ = profile(N)
    final = march!(f₀, scheme, CMP_C, round(Int, N/CMP_C))
    return sqrt(Δx*sum(abs2, final .- f₀))
end

# `fmax = 2.5` brackets every profile below: the sine peaks at 1.5, the gaussian
# at 2.0 and the pulse at 1.5.
cmp_schemes() = [
    ("Upwind",                   Upwind()),
    ("Godunov constant",         Godunov(PiecewiseConstant())),
    ("SemiLagrangian linear",    SemiLagrangian(LinearSpline())),
    ("LaxWendroff",              LaxWendroff()),
    ("Godunov VanLeer",          Godunov(PiecewiseLinear(), VanLeer())),
    ("SemiLagrangian quadratic", SemiLagrangian(QuadraticSpline())),
    ("SemiLagrangian cubic",     SemiLagrangian(CubicSpline())),
    ("PFC",                      PFC(fmin = 0.0, fmax = 2.5)),
]

@testset "Scheme comparison" begin

    err = Dict{Tuple{String,Symbol},Float64}()
    for (name, scheme) in cmp_schemes()
        err[(name, :sine)]   = traversal_error(scheme, cmp_sine)
        err[(name, :gauss)]  = traversal_error(scheme, cmp_gauss)
        err[(name, :square)] = traversal_error(scheme, cmp_square)
    end

    println("  L² error after one traversal, N = ", CMP_N, ", c = ", CMP_C)
    println("  ", rpad("scheme", 26), lpad("sine", 12), lpad("gaussian", 12), lpad("square", 12))
    for (name, _) in cmp_schemes()
        println("  ", rpad(name, 26),
                lpad(round(err[(name, :sine)];   sigdigits = 4), 12),
                lpad(round(err[(name, :gauss)];  sigdigits = 4), 12),
                lpad(round(err[(name, :square)]; sigdigits = 4), 12))
    end

    @testset "the equivalent schemes agree in the error too" begin
        # `test_amplification` shows these three share an amplification factor
        # mode by mode. That is the reason; this is the consequence, after 1280
        # steps rather than one. Measured max|Δ| in the final state: 7.3e-15
        # against Godunov constant, 2.9e-14 against the linear spline -- not
        # bit-identical, because the spline reaches the same value by a
        # different route, but far below any error being compared here.
        for profile in (:sine, :gauss, :square)
            @test err[("Godunov constant", profile)] ≈ err[("Upwind", profile)] rtol = 1e-9
            @test err[("SemiLagrangian linear", profile)] ≈ err[("Upwind", profile)] rtol = 1e-9
        end
    end

    @testset "on smooth data, order is what matters" begin
        # Strictly ordered by the order of the scheme, over five decades:
        #
        #   SemiLagrangian cubic      2.46e-08
        #   PFC                       2.30e-07
        #   SemiLagrangian quadratic  5.02e-06
        #   LaxWendroff               4.68e-05
        #   Godunov VanLeer           7.03e-03
        #   Upwind                    8.08e-03
        #
        # `Godunov VanLeer` sits second from the bottom, barely ahead of upwind,
        # because the limiter clips smooth extrema -- `test_convergence` holds it
        # to first order in L¹ for exactly this reason. Remember where it is:
        # the next testset is the same list on a discontinuity.
        smooth = [err[(n, :sine)] for n in
                  ("SemiLagrangian cubic", "PFC", "SemiLagrangian quadratic",
                   "LaxWendroff", "Godunov VanLeer", "Upwind")]
        @test issorted(smooth)
        # and the spread is five decades, not a tie being reported as an order
        @test smooth[end]/smooth[1] > 1e5
    end

    @testset "on a discontinuity, the ranking inverts" begin
        # Godunov's theorem, as a measurement. On the square pulse:
        #
        #   Godunov VanLeer           8.45e-03    <- best
        #   SemiLagrangian cubic      2.01e-02
        #   PFC                       2.65e-02
        #   SemiLagrangian quadratic  2.72e-02
        #   LaxWendroff               4.53e-02
        #   Upwind                    6.32e-02
        #
        # The scheme that was second from the bottom on the sine is now first,
        # and the one that led by five decades is second. Nothing else in the
        # suite says this, and it is the single most useful thing to know when
        # choosing a scheme: there is no ordering of these that survives a
        # change of problem class.
        @test err[("Godunov VanLeer", :square)] < err[("SemiLagrangian cubic", :square)]
        @test err[("Godunov VanLeer", :square)] < err[("PFC", :square)]
        @test err[("Godunov VanLeer", :square)] < err[("LaxWendroff", :square)]

        # The inversion, quantified: a ratio of 2.9e5 one way becomes 0.42 the
        # other, which is a swing of six decades in relative standing.
        smooth_ratio = err[("Godunov VanLeer", :sine)]/err[("SemiLagrangian cubic", :sine)]
        rough_ratio  = err[("Godunov VanLeer", :square)]/err[("SemiLagrangian cubic", :square)]
        println("  Godunov VanLeer / SemiLagrangian cubic: ",
                round(smooth_ratio; sigdigits = 3), " on the sine, ",
                round(rough_ratio; sigdigits = 3), " on the pulse")
        @test smooth_ratio > 1e4
        @test rough_ratio < 1.0
    end

    @testset "PFC beats the quadratic spline where it is not a tie" begin
        # More accurate on both smooth profiles *and* an order of magnitude
        # cheaper -- 3.5 ns per cell per step against 65 to 85, the spline
        # prefilter allocating where PFC does not (`test_allocations` pins
        # both). The work-precision report has the quadratic spline dominated on
        # all three profiles once cost is counted.
        #
        # Only the two smooth ratios are asserted. On the square pulse the two
        # are within 2.6% of each other at this resolution -- 2.646e-2 against
        # 2.717e-2 -- and the sign of a 2.6% gap is not a property worth gating:
        # it is already the other way round at N = 256, where PFC reads 3.416e-2
        # against the spline's 3.392e-2. Asserting it would buy a flapping test
        # and no information.
        @test err[("PFC", :sine)]  < err[("SemiLagrangian quadratic", :sine)]/10
        @test err[("PFC", :gauss)] < err[("SemiLagrangian quadratic", :gauss)]/2
        println("  PFC / SemiLagrangian quadratic: ",
                join([string(profile, " ",
                             round(err[("PFC", profile)]/err[("SemiLagrangian quadratic", profile)];
                                   sigdigits = 3))
                      for profile in (:sine, :gauss, :square)],
                      ", "), "   (the pulse is a tie, and not asserted)")
    end
end

# ---------------------------------------------------------------- complexity

"""
    scaling_exponent(work, sizes; budget = 0.05, minreps = 5)

Least-squares slope of `log(time)` against `log(N)`: the `p` in `t ∝ Nᵖ`.

Each size is warmed twice before it is clocked and the **minimum** over the
samples is taken. Both matter. Julia specialises the kernel per operator type,
so a first call times the compiler -- a mistake this repository has already made
once in `verification/scheme-comparison.jl`, where it reported upwind as slower
than `PFCNonUniform` when it is about twice as fast. And the minimum, rather
than the mean, is the estimator that a shared machine's interruptions cannot
inflate.

Sampling is by **time budget per size, not by a fixed count**. The two are not
the same thing where it matters. A fitted slope is most sensitive to its
endpoints, and the smallest size is both the shortest measurement -- the one a
scheduler interruption can inflate by a multiple rather than a percent -- and
the one where inflation flattens the fit toward zero. A fixed count spends the
same number of samples on the point that needs many and the point that needs
few; here the cheap points get thousands and the expensive ones get `minreps`,
which is the allocation the estimator actually wants. It also holds the cost of
this testset flat at roughly `budget` per size regardless of the machine and of
whether the run carries `--check-bounds=yes` or coverage instrumentation, both
of which move these kernels by an order of magnitude.
"""
function scaling_exponent(work, sizes; budget = 0.05, minreps = 5)
    times = map(sizes) do n
        f = work(n)
        f(); f()                                   # compile, off the clock
        t = Inf
        spent = 0.0
        reps = 0
        while reps < minreps || spent < budget
            t0 = time_ns()
            f()
            δ = (time_ns() - t0)/1e9
            t = min(t, δ)
            spent += δ
            reps += 1
        end
        t
    end
    X = log.(collect(float.(sizes)))
    y = log.(times)
    n = length(X)
    p = (n*sum(X.*y) - sum(X)*sum(y))/(n*sum(X.^2) - sum(X)^2)
    return p, times
end

@testset "Collision operators: complexity class" begin
    # The one timing statement in the suite that is a gate rather than a report.
    #
    # It is gateable because the quantity is a complexity class, not a duration:
    # `BGK` doubles when N doubles and `Landau1P` quadruples, so the exponents
    # are 1 and 2 with nothing in between for noise to land on.
    #
    # **Measured in the mode this actually runs in.** `Pkg.test()` passes
    # `--check-bounds=yes`, so these kernels are about three times slower here
    # than the same code timed from a plain `julia --project=.`, and quoting the
    # faster figures would leave anyone debugging a failure comparing against
    # numbers the suite never prints. Five trials under `--check-bounds=yes`:
    #
    #   BGK       13.6 27.4 55.0 110 us    ratios 2.00-2.09   exponent 1.003-1.028
    #   Landau1P  102 401 1600 6385 us     ratios 3.91-4.07   exponent 1.990-2.001
    #
    # Held to a threshold at 1.5, which is some forty measurement errors from
    # either. The `Coverage` job runs the same assertions under
    # `--code-coverage=user`, where every number above grows by a factor of
    # about fourteen: 194 to 1551 us and 673 us to 43 ms. The exponents come out
    # 0.983-1.005 and 1.997-2.008 -- unmoved, because a uniform slowdown cancels
    # in a ratio, which is the whole reason this quantity can be gated where a
    # duration cannot.
    #
    # That is also why the 10 us floor in `runbenchmarks.jl` does not apply. The
    # floor is there because comparing a *duration* against a baseline stored on
    # another day needs a tolerance that survives 20-50% run-to-run drift, and
    # below 10 us it does not -- that file measures a 10 us floor still tripping
    # once in four runs. Nothing here is compared against a stored number: the
    # four measurements are taken seconds apart in one process, and only their
    # ratio is read.
    #
    # It catches a real regression that nothing else would: an accidental O(N²)
    # in `BGK` -- a moment recomputed inside the velocity loop, say -- leaves the
    # allocation gate happy and every physics assertion passing, and shows up
    # only as a machine that has become slow.
    #
    # **The advection kernels are deliberately not gated this way**, and that is
    # a measurement rather than an omission: their per-call times are
    # sub-microsecond, where timer resolution and cache state dominate, and a
    # fourfold size increase on an O(N) kernel measured ratios from 3.29 to 8.00
    # against the 4 it should give. `benchmark/workprecision.jl` reports their
    # cost without asserting it.
    function collision_work(op)
        return function (n)
            v = collect(range(-4, 4; length = n))
            src = @. exp(-v^2)
            dst = similar(src)
            ws = collision_workspace(op, n)
            return () -> collide!(dst, src, op, v, 0.1, ws)
        end
    end

    p_bgk, t_bgk = scaling_exponent(collision_work(BGK(1e-2)),
                                    (800, 1600, 3200, 6400))
    p_landau, t_landau = scaling_exponent(collision_work(Landau1P(1e-2)),
                                          (100, 200, 400, 800))

    println("  BGK       N = 800..6400  ",
            join([string(round(t*1e6; digits = 1), "us") for t in t_bgk], " "),
            "   exponent ", round(p_bgk; digits = 2))
    println("  Landau1P  N = 100..800   ",
            join([string(round(t*1e6; digits = 1), "us") for t in t_landau], " "),
            "   exponent ", round(p_landau; digits = 2))

    @test p_bgk < 1.5           # linear in the velocity grid
    @test p_landau > 1.5        # the double loop over velocity pairs
    @test p_landau - p_bgk > 0.5
end
