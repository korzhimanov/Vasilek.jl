using Vasilek

include(joinpath(@__DIR__, "scheme_cases.jl"))

"""
Direction symmetry, and the degenerate cases either side of it.

`test_golden`, `test_convergence` and `test_invariants` all run at c = 0.4 > 0,
and `c > 0` and `c < 0` are separate branches in every scheme that has one.
Only `test_1d_advection` touched the negative branch at all, and only for a
single step against a hand-chosen tolerance.

Rather than duplicate four suites at negative c, this asserts the identity that
ties the two branches together. Reversing the grid reverses the flow direction,
so a scheme whose branches agree must commute with the reversal:

    R∘A(+c) == A(−c)∘R,        R(f)[i] = f[n+1−i]

That this is enough was checked rather than assumed: the backward-direction
convergence orders and invariants were measured and agree with the forward ones
to three decimals, which is exactly what the identity predicts. An index swapped
by one in either branch breaks it immediately.

**Why some schemes are exact and others are not.** For `Upwind` and `Godunov`
the mirrored expression is the *same* floating-point expression on the same
operands in the same order, so the identity holds bit-for-bit. `LaxWendroff`
evaluates `f[i⁺] - 2f[i] + f[i⁻]` left to right, and mirroring exchanges `f[i⁺]`
with `f[i⁻]`, which reassociates the sum — hence two ULP. The semi-Lagrangian
family goes through a global periodic prefilter, where a ULP at one end reaches
the other, hence 4.5e-14.

Measured over five datasets × four Courant numbers.
"""
const SYM_N = 128

"Structurally varied, and deterministic: no RNG, so no new test dependency."
sym_data() = [
    ("smooth",  [1.0 + 0.5*sin(2π*(i-1)/SYM_N) for i = 1:SYM_N]),
    ("pulse",   [0.4 < (i-1)/SYM_N < 0.6 ? 1.0 : 0.0 for i = 1:SYM_N]),
    ("kink",    [abs((i-1)/SYM_N - 0.5) for i = 1:SYM_N]),
    ("rough",   [sin(i^2/7.0)*cos(3.0*i) for i = 1:SYM_N]),
    ("nyquist", [0.1 + (iseven(i) ? 1.0 : -1.0) for i = 1:SYM_N]),
]

"""
The canonical scheme list, with `PFC` bracketed to the data it will be given.

`PFC` requires global initial bounds by contract, and these datasets have
different ranges, so the bounds have to come from the data rather than from a
constant shared across the suite.
"""
function sym_schemes(f)
    lo, hi = minimum(f), maximum(f)
    pad = 1e-9*max(1.0, abs(lo), abs(hi))
    return uniform_schemes(fmin = lo - pad, fmax = hi + pad)
end

"Schemes whose mirrored expression is operand-for-operand identical."
const SYM_EXACT = ("Upwind", "Godunov_constant", "Godunov_linear",
                   "Godunov_linear_VanLeer")

sym_tol(name) = name in SYM_EXACT           ? 0.0   :   # measured 0.0
                startswith(name, "SemiLag") ? 1e-12 :   # measured 4.5e-14
                                              1e-14     # measured 4.4e-16

step1(f, scheme, c) = march!(f, scheme, c, 1)

@testset "Direction symmetry: R∘A(+c) == A(−c)∘R" begin
    worst = Dict{String,Float64}()
    for (_, f) in sym_data(), c in (0.4, 0.7, 0.15, 0.9)
        for (name, scheme) in sym_schemes(f)
            a = reverse(step1(f, scheme, c))
            b = step1(reverse(f), scheme, -c)
            d = maximum(abs, a .- b)
            worst[name] = max(get(worst, name, 0.0), d)
            @test d ≤ sym_tol(name)
        end
    end
    for (name, _) in sym_schemes(sym_data()[1][2])
        println("  mirror ", rpad(name, 26), worst[name],
                "   (tol ", sym_tol(name), ")")
    end
end

@testset "A constant field stays constant" begin
    # Consistency: every scheme's stencil weights sum to one, so a uniform
    # field is a fixed point. This is the cheapest test that catches a
    # mis-weighted flux, and it exercises `Godunov._ratio` on the exactly-zero
    # difference its `≈ 0.0` branch is written for.
    #
    # Measured worst relative deviation: 0.0 for Upwind, LaxWendroff and the
    # linear spline; 1.9e-16 for Godunov and PFC; 4.4e-16 for the cubic spline.
    for val in (0.3, 1.0, 7.5), c in (0.4, -0.4, 0.9, -0.15)
        f = fill(val, SYM_N)
        for (name, scheme) in sym_schemes(f)
            out = step1(f, scheme, c)
            @test maximum(abs, out .- val)/val ≤ 1e-14
        end
    end
end

@testset "c = 0 is the identity" begin
    # Note this runs the *negative* branch of every scheme that has one:
    # `c > 0` is false at zero. Exact everywhere except the quadratic and cubic
    # splines, where the prefilter round-trip costs 4.4e-16 and 7.8e-16.
    for (label, f) in sym_data()
        for (name, scheme) in sym_schemes(f)
            out = step1(f, scheme, 0.0)
            tol = startswith(name, "SemiLagrangian_qu") ||
                  startswith(name, "SemiLagrangian_cu") ? 1e-14 : 0.0
            @test maximum(abs, out .- f) ≤ tol
        end
    end
end

@testset "SemiLagrangian is periodic in the Courant number" begin
    # Backward characteristic tracing on a periodic grid cannot distinguish a
    # displacement of c from one of c ± N, so the two must agree. This is the
    # only test of the `extrapolate(itp, Periodic(OnCell()))` path in
    # `_sample!`, which the interior branch of the loop never reaches.
    #
    # Measured worst: 6.4e-14.
    splines = [("linear", LinearSpline()), ("quadratic", QuadraticSpline()),
               ("cubic", CubicSpline())]
    worst = 0.0
    for (_, f) in sym_data(), c in (0.4, -0.4, 0.85, -0.85)
        for (_, sp) in splines
            s = SemiLagrangian(sp)
            a = step1(f, s, c)
            worst = max(worst, maximum(abs, a .- step1(f, s, c + SYM_N)),
                               maximum(abs, a .- step1(f, s, c - SYM_N)))
            @test a ≈ step1(f, s, c + SYM_N) atol=1e-12
            @test a ≈ step1(f, s, c - SYM_N) atol=1e-12
        end
    end
    println("  SemiLagrangian A(c) vs A(c ± N): worst max|Δ| = ", worst)
end

@testset "SemiLagrangian is stable beyond the Courant limit" begin
    # Semi-Lagrangian tracing has no CFL condition -- that is its whole point,
    # and nothing asserted it. Every other scheme here blows up above |c| = 1.
    #
    # Measured max|f| after 500 steps at c = 3.7, from an initial 1.5:
    # 1.4406 (linear, diffusive), 1.49999 (quadratic), 1.49999 (cubic).
    f₀ = [1.0 + 0.5*sin(2π*(i-1)/SYM_N) for i = 1:SYM_N]
    for (label, sp) in [("linear", LinearSpline()), ("quadratic", QuadraticSpline()),
                        ("cubic", CubicSpline())]
        f = march!(f₀, SemiLagrangian(sp), 3.7, 500)
        println("  SemiLagrangian ", rpad(label, 10), " at c = 3.7, 500 steps: max|f| = ",
                maximum(abs, f))
        @test maximum(abs, f) ≤ maximum(abs, f₀) + 1e-12
        @test all(isfinite, f)
    end
end
