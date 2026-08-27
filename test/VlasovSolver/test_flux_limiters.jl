using Vasilek

@testset "Flux limiters" begin
    # Limiters are callable singleton types now, so an unrecognised one is a
    # MethodError where it is written rather than a call on `nothing` from
    # inside the hot loop.
    for limiter in (VanLeer(),)
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

    # NoLimiter is the identity that makes PiecewiseLinear unstable
    @test NoLimiter()(0.0) == 1.0
    @test NoLimiter()(17.0) == 1.0

    # options that do not exist fail where they are written
    @test_throws MethodError Godunov(:Riemann_linear)
    @test_throws MethodError Godunov(PiecewiseLinear(), :VanLeer)
end
