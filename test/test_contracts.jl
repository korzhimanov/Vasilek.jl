using Vasilek

include(joinpath(@__DIR__, "scheme_cases.jl"))

# API contract tests: the promises the types make outside the numerics.
#
# `PFC`'s bounds carry a whole paragraph of docstring and the memory of a bug
# that made two overloads disagree by a factor of 740, and `checked` exists
# solely to police them -- yet nothing exercised either the check or the type
# parameter that compiles it away.
@testset "PFC bounds check" begin
    N = 64
    f = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]      # ranges over [0.5, 1.5]

    @testset "the check fires, in both branches" begin
        # Both signs, because the assertion sits above the branch and a future
        # edit could easily move it inside one of them.
        over = copy(f);  over[10] = 5.0
        under = copy(f); under[10] = -1.0
        for c in (0.4, -0.4)
            @test_throws AssertionError march!(over, PFC(fmin = 0.0, fmax = 2.0), c, 1)
            @test_throws AssertionError march!(under, PFC(fmin = 0.0, fmax = 2.0), c, 1)
        end

        # Touching the bounds exactly is legal: the assertions are `≤`.
        exact = copy(f); exact[1] = 0.0; exact[2] = 2.0
        @test march!(exact, PFC(fmin = 0.0, fmax = 2.0), 0.4, 1) isa Vector{Float64}
    end

    @testset "checked = false compiles it away without changing the answer" begin
        # The point of the type parameter: identical numerics on valid data,
        # no minimum/maximum pass. Bit-for-bit, not approximately.
        for c in (0.4, -0.4)
            on  = march!(f, PFC(fmin = 0.0, fmax = 2.0, checked = true), c, 4)
            off = march!(f, PFC(fmin = 0.0, fmax = 2.0, checked = false), c, 4)
            @test on == off
        end

        # And on data outside the bounds it does not throw -- it produces a
        # wrong answer quietly, which is precisely the trade the flag makes.
        over = copy(f); over[10] = 5.0
        @test march!(over, PFC(fmin = 0.0, fmax = 2.0, checked = false), 0.4, 1) isa Vector{Float64}
    end

    @testset "the flag and the element type are type parameters" begin
        # `checked` has to be a type parameter for the branch to be elided, and
        # the element type has to survive `promote`, or a Float32 run silently
        # widens.
        @test PFC(fmin = 0.0, fmax = 2.0)                    isa PFC{Float64, true}
        @test PFC(fmin = 0.0, fmax = 2.0, checked = false)   isa PFC{Float64, false}
        @test PFC(fmin = 0, fmax = 2)                        isa PFC{Float64, true}
        @test PFC(fmin = 0.0f0, fmax = 2.0f0)                isa PFC{Float32, true}
        # `float(::Int)` is Float64, so integer bounds *widen* a Float32
        # pairing rather than narrowing to it. Pinned because it is the
        # surprising way round.
        @test PFC(fmin = 0, fmax = 2.0f0)                    isa PFC{Float64, true}
    end
end
