using Vasilek
using LinearAlgebra
# Imported, not `using`: Interpolations and Vasilek both export `OnGrid`, and a
# test run puts every file in `Main`.
import Interpolations

"""
The semi-Lagrangian splines are the package's own: a periodic prefilter solved
by Sherman–Morrison against a factorisation the workspace holds, and an
evaluation that follows Interpolations.jl 0.16, which `SemiLagrangian` called
until 0.2. This holds the two to each other.

  * The linear spline is the same to the bit, every cell, the periodic map
    included.
  * The quadratic and cubic splines evaluate with Interpolations' own weights
    in its own order, so the coefficients alone differ: the same system, solved
    by a different factorisation. Measured worst, relative to the data's
    maximum, over the cases below: 4.4e-15 for both.

What the prefilter has to be is checked on its own as well: the spline it
builds passes through the data at the knots, and its coefficients solve the
periodic system.
"""

# The 0.2 kernel: Interpolations' periodic prefilter and extrapolation, called
# per point as it was. `interpolate!` prefilters its argument in place, hence a
# fresh buffer per call, which for the linear spline repeats the first point.
spline_itp(::LinearSpline) = Interpolations.BSpline(Interpolations.Linear())
spline_itp(::QuadraticSpline) =
    Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Periodic(Interpolations.OnCell())))
spline_itp(::CubicSpline) =
    Interpolations.BSpline(Interpolations.Cubic(Interpolations.Periodic(Interpolations.OnCell())))

function spline_reference(src, spline, c)
    buf = spline isa LinearSpline ? vcat(src, src[1]) : copy(src)
    itp = Interpolations.interpolate!(buf, spline_itp(spline))
    etp = Interpolations.extrapolate(itp, Interpolations.Periodic(Interpolations.OnCell()))
    return [1 ≤ i - c ≤ length(itp) ? itp(i - c) : etp(i - c) for i in eachindex(src)]
end

"`test_symmetry.jl`'s five data sets, at any `n`."
spline_data(n) = [
    ("smooth",  [1.0 + 0.5*sin(2π*(i-1)/n) for i = 1:n]),
    ("pulse",   [0.4 < (i-1)/n < 0.6 ? 1.0 : 0.0 for i = 1:n]),
    ("kink",    [abs((i-1)/n - 0.5) for i = 1:n]),
    ("rough",   [sin(i^2/7.0)*cos(3.0*i) for i = 1:n]),
    ("nyquist", [0.1 + (iseven(i) ? 1.0 : -1.0) for i = 1:n]),
]

# Inside and beyond one cell, both ways, and past the Courant limit; then the
# seams. 0, ±1 and 2 put points on knots, the linear spline's last one among
# them; ±(n + 0.4) sends every point through the periodic map; and the last two
# make `mod` round up to the period, which puts a point on the upper end of the
# knots: `n + 1` for the linear spline (`eps()/2`), `n + 1/2` for the others
# (`nextfloat(0.5)`).
spline_courants(n) = (0.4, -0.77, 1.2, 3.7, 0.15,
                      0.0, 1.0, -1.0, 2.0, n + 0.4, -(n + 0.4), eps()/2, nextfloat(0.5))

const SPLINES = (("linear", LinearSpline()), ("quadratic", QuadraticSpline()),
                 ("cubic", CubicSpline()))

@testset "SemiLagrangian splines" begin
    @testset "against Interpolations" begin
        worst = Dict("quadratic" => 0.0, "cubic" => 0.0)
        for n in (7, 64, 161), (_, f) in spline_data(n), c in spline_courants(n),
            (label, sp) in SPLINES
            s = SemiLagrangian(sp)
            out = advect!(similar(f), f, s, c, workspace(s, n))
            ref = spline_reference(f, sp, c)
            if sp isa LinearSpline
                @test out == ref
            else
                d = maximum(abs, out .- ref)/maximum(abs, ref)
                worst[label] = max(worst[label], d)
                @test d ≤ 1e-13
            end
        end
        println("  SemiLagrangian vs Interpolations: worst max|Δ|/max|f| = ", worst)
    end

    @testset "the interior and the edges are evaluated alike" begin
        # `advect!` takes the cells whose stencil lies inside the knots through
        # a loop with no periodic map and no wrapped index; the general loop,
        # over every cell, has to give the same bits.
        for n in (7, 64), (_, f) in spline_data(n), c in spline_courants(n), (_, sp) in SPLINES
            s = SemiLagrangian(sp)
            ws = workspace(s, n)
            out = advect!(similar(f), f, s, c, ws)
            @test out == Vasilek.Advection._sample!(similar(f), ws.coefficients, sp,
                                                    float(c), n, 1:n)
        end
    end

    @testset "the prefilter interpolates, and solves the periodic system" begin
        # From three cells, where every pair of cells neighbours each other.
        for n in (3, 4, 7, 64), (_, f) in spline_data(n),
            sp in (QuadraticSpline(), CubicSpline())
            scale = maximum(abs, f)
            ws = workspace(SemiLagrangian(sp), n)
            coefficients = Vasilek.Advection._prefilter!(ws, sp, f)
            at_knots = [Vasilek.Advection._evaluate(coefficients, n, sp, Float64(i), Val(true))
                        for i = 1:n]
            @test maximum(abs, at_knots .- f) ≤ 1e-14*scale
            a, b = Vasilek.Advection._diagonals(Float64, sp)
            M = diagm(0 => fill(a, n), 1 => fill(b, n - 1), -1 => fill(b, n - 1))
            M[1, n] += b
            M[n, 1] += b
            @test maximum(abs, M*coefficients .- f) ≤ 1e-14*scale
        end
    end

    @testset "a workspace belongs to its spline" begin
        # Another spline's workspace has the right length, and holds another
        # matrix's factorisation: it would give a plausible wrong answer.
        n = 64
        f = spline_data(n)[1][2]
        for (_, sp) in SPLINES, (_, other) in SPLINES
            sp == other && continue
            s = SemiLagrangian(sp)
            @test_throws ArgumentError advect!(similar(f), f, s, 0.4, workspace(SemiLagrangian(other), n))
        end
        s = SemiLagrangian(CubicSpline())
        @test workspace(s, n, Float32) isa Vasilek.Advection.SplineWorkspace{Float32, CubicSpline}
        # The factorisation is read on every step; a dest on top of it would
        # change every later one.
        ws = workspace(s, n)
        @test_throws ArgumentError advect!(view(ws.z, 1:n), f, s, 0.4, ws)
        @test_throws ArgumentError advect!(view(ws.rdiag, 1:n), f, s, 0.4, ws)
    end

    @testset "the quadratic spline in single precision" begin
        # `test_contracts.jl` holds the cubic and the linear spline to this.
        n = 64
        g = Float32[1.0 + 0.5*sin(2π*(i-1)/n) for i = 1:n]
        s = SemiLagrangian(QuadraticSpline())
        dst = similar(g)
        @test (@inferred advect!(dst, g, s, 0.4f0, workspace(s, n, Float32))) === dst
        @test maximum(abs, dst .- advect!(similar(g, Float64), Float64.(g), s, 0.4)) < 1e-5
    end

    @testset "a step computes in the data's type" begin
        # Whatever the type of `c`: it used to set the type of every point and
        # weight. A Float16 `c` rounded `i - c` past the interior's last cell
        # (n = 5000) and the period past the last knot (n = 2051), and the
        # stencil read beyond the coefficients under @inbounds; a Float64 `c`
        # left BigFloat data with Float64-accurate weights.
        for (_, sp) in SPLINES
            s = SemiLagrangian(sp)
            for n in (2051, 5000)
                f = [1 + 0.5sin(2π*i/n) + 0.1cos(6π*i/n) for i = 1:n]
                c = Float16(0.4)
                @test advect!(similar(f), f, s, c) == advect!(similar(f), f, s, Float64(c))
            end
            g = Float32[1 + 0.5sin(2π*i/64) for i = 1:64]
            @test advect!(similar(g), g, s, 0.4) == advect!(similar(g), g, s, 0.4f0)
            gb = big.(g)
            # 2.7e-17 for the cubic spline, the gap between 0.4 and big"0.4";
            # 1.7e-15 when the weights were Float64
            @test maximum(abs, advect!(similar(gb), gb, s, 0.4) .-
                               advect!(similar(gb), gb, s, big"0.4")) < 1e-16
        end
    end

    @testset "a type that cannot hold the step is refused" begin
        s = SemiLagrangian(CubicSpline())
        # Float16 holds the integers only up to 2048, and the points and the
        # period need every knot exact; a wider workspace takes the data.
        h = Float16[1 + 0.5sin(2π*i/5000) for i = 1:5000]
        @test_throws ArgumentError advect!(similar(h), h, s, 0.4)
        @test all(isfinite, advect!(similar(h), h, s, 0.4, workspace(s, 5000, Float32)))
        h = Float16[1 + 0.5sin(2π*i/64) for i = 1:64]
        @test_throws DomainError advect!(similar(h), h, s, 1e6)
        # A narrower workspace would round the data on the way in.
        f = [1 + 0.5sin(2π*i/64) for i = 1:64]
        @test_throws ArgumentError advect!(similar(f), f, s, 0.4, workspace(s, 64, Float32))
    end

    @testset "the workspace is the spline's, in a float type" begin
        # A scheme with an abstract parameter labelled its workspace with it
        # and then refused it; an Int workspace threw on its reciprocal pivots.
        f = [1 + 0.5sin(2π*i/64) for i = 1:64]
        for (_, sp) in SPLINES
            s = SemiLagrangian{Vasilek.Advection.AbstractSpline}(sp)
            @test advect!(similar(f), f, s, 0.4, workspace(s, 64)) ==
                  advect!(similar(f), f, SemiLagrangian(sp), 0.4)
            ws = workspace(SemiLagrangian(sp), 64, Int)
            @test ws isa Vasilek.Advection.SplineWorkspace{Float64, typeof(sp)}
            k = collect(1:64)
            @test advect!(similar(f), k, SemiLagrangian(sp), 0.4, ws) ==
                  advect!(similar(f), float.(k), SemiLagrangian(sp), 0.4)
        end
    end
end
