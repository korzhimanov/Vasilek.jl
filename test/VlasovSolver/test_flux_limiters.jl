const LIMITERS = (:VanLeer,)

include(joinpath(dirname(@__FILE__),"..","..","src","VlasovSolver","Limiters.jl"))

function test_limiter(limiter_name)
    limiter = Limiters.generate_limiter(limiter_name)

    # test special points
    @test limiter(0.0) ≈ 0.0
    @test limiter(1.0) ≈ 1.0

    # test symmetricity
    for r in (2.0, 3.0, 5.0, 11.0)
        @test limiter(r) ≈ r*limiter(1/r)
    end

    # test TVD condition
    for r in range(0.1, 0.9; step=0.1)
        @test r ≤ limiter(r) ≤ 2*r
    end
    for r in range(1.1, 1.9; step=0.1)
        @test 1 ≤ limiter(r) ≤ r
    end
    for r in range(2.0, 10.0; step=1.0)
        @test 1.0 ≤ limiter(r) ≤ 2.0
    end
end

for limiter_name in LIMITERS
    test_limiter(limiter_name)
end
