using Vasilek

@testset "Flux limiters" begin
    # Limiters are callable singleton types now, so an unrecognised one is a
    # MethodError where it is written rather than a call on `nothing` from
    # inside the hot loop.
    for limiter in (VanLeer(), Superbee())
        # special points
        @test limiter(0.0) ≈ 0.0
        @test limiter(1.0) ≈ 1.0

        # symmetry
        for r in (2.0, 3.0, 5.0, 11.0)
            @test limiter(r) ≈ r*limiter(1/r)
        end

        # TVD region
        for r in range(0.1, 0.9; step = 0.1)
            @test r ≤ limiter(r) ≤ 2*r
        end
        for r in range(1.1, 1.9; step = 0.1)
            @test 1 ≤ limiter(r) ≤ r
        end
        for r in range(2.0, 10.0; step = 1.0)
            @test 1.0 ≤ limiter(r) ≤ 2.0
        end
    end

    # Superbee is the upper edge of that region, piece by piece: 2r, then 1,
    # then r, then 2. Every limiter in the region lies at or below it, VanLeer
    # included, and both are zero where the slopes change sign.
    sb = Superbee()
    for r in range(0.0, 0.5; length = 11)
        @test sb(r) == 2r
    end
    for r in range(0.5, 1.0; length = 11)
        @test sb(r) == 1.0
    end
    for r in range(1.0, 2.0; length = 11)
        @test sb(r) == r
    end
    for r in (2.0, 3.5, 10.0, 1e6)
        @test sb(r) == 2.0
    end
    for r in range(0.0, 10.0; length = 101)
        @test VanLeer()(r) ≤ sb(r)
    end
    for r in (-1e6, -1.0, -0.3, -0.0)
        @test sb(r) == 0.0
        @test VanLeer()(r) == 0.0
    end

    # NoLimiter is the identity, with which Godunov(PiecewiseLinear()) is
    # LaxWendroff (`test_convergence`, `test_amplification`)
    @test NoLimiter()(0.0) == 1.0
    @test NoLimiter()(17.0) == 1.0

    # options that do not exist fail where they are written
    @test_throws MethodError Godunov(:Riemann_linear)
    @test_throws MethodError Godunov(PiecewiseLinear(), :VanLeer)
end
