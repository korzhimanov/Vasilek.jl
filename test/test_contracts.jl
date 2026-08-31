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

# The argument checks `advect!` runs, and the three ways of calling it wrongly
# that used to produce a plausible answer instead of an error. Each case below
# was measured misbehaving before the check existed; the numbers are in the
# `_validate` docstring.
@testset "advect! argument validation" begin
    N = 64
    f = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
    all_schemes = vcat(uniform_schemes(fmin = 0.0, fmax = 2.0),
                       [("PFCNonUniform", PFCNonUniform(fill(0.05, N); fmin = 0.0, fmax = 2.0))])

    @testset "dest === src is rejected" begin
        # Silently wrong before: 0.013 for Upwind, 0.21 for PFC. SemiLagrangian
        # and PFCNonUniform happened to survive it, by holding a scratch copy --
        # an accident, not a contract, so all schemes reject it alike.
        for (name, scheme) in all_schemes
            a = copy(f)
            c = name == "PFCNonUniform" ? 0.02 : 0.4
            @test_throws ArgumentError advect!(a, a, scheme, c, workspace(scheme, N))
        end
    end

    @testset "unequal lengths are rejected" begin
        # Was half-loud: length(dest) < length(src) truncated silently in the
        # finite-difference schemes, wrapping periodically at the shorter
        # length, while the other direction threw BoundsError.
        for (name, scheme) in uniform_schemes(fmin = 0.0, fmax = 2.0)
            @test_throws DimensionMismatch advect!(Vector{Float64}(undef, 32), f,
                                                   scheme, 0.4, workspace(scheme, N))
            @test_throws DimensionMismatch advect!(Vector{Float64}(undef, N), f[1:32],
                                                   scheme, 0.4, workspace(scheme, N))
        end
    end

    @testset "a workspace of the wrong size is rejected" begin
        # Oversized was the dangerous one: `interpolate!` prefilters the whole
        # buffer, so a SemiLagrangian handed a longer workspace returned
        # garbage (max|Δ| ≈ 0.4) without complaint. Undersized already threw
        # BoundsError, which is loud but from the wrong place.
        for spline in (LinearSpline(), QuadraticSpline(), CubicSpline())
            s = SemiLagrangian(spline)
            dst = similar(f)
            @test_throws DimensionMismatch advect!(dst, f, s, 0.4, workspace(s, N + 10))
            @test_throws DimensionMismatch advect!(dst, f, s, 0.4, workspace(s, N - 10))
            @test advect!(dst, f, s, 0.4, workspace(s, N)) === dst    # the right one works
        end

        s = PFCNonUniform(fill(0.05, N); fmin = 0.0, fmax = 2.0)
        dst = similar(f)
        @test_throws DimensionMismatch advect!(dst, f, s, 0.02, workspace(s, N + 10))
        @test_throws DimensionMismatch advect!(dst, f, s, 0.02, workspace(s, N - 10))

        # and a scheme whose grid does not match the data
        wrong = PFCNonUniform(fill(0.05, 32); fmin = 0.0, fmax = 2.0)
        @test_throws DimensionMismatch advect!(dst, f, wrong, 0.02, workspace(wrong, 32))
    end

    @testset "the minimum problem size is uniform" begin
        # `Godunov(PiecewiseLinear())` reaches two neighbours either side, so it
        # threw BoundsError at n = 2 while every other scheme quietly produced a
        # degenerate answer. Now all of them agree on where the floor is.
        for (name, scheme) in uniform_schemes(fmin = -1.0, fmax = 2.0)
            @test_throws ArgumentError march!([1.0, 1.5], scheme, 0.4, 1)
        end
        for n = 3:8
            g = [1.0 + 0.3*sin(2π*(i-1)/n) for i = 1:n]
            for (name, scheme) in uniform_schemes(fmin = -1.0, fmax = 2.0), c in (0.4, -0.4)
                out = march!(g, scheme, c, 1)
                @test length(out) == n
                @test abs(sum(out) - sum(g))/abs(sum(g)) < 1e-13
            end
        end
    end
end

@testset "a workspace carries no state between calls" begin
    # `workspace` hands back `undef` memory, so a scheme that ever read a slot
    # before writing it would give a different answer on a reused workspace
    # than on a fresh one -- and would do so nondeterministically. Nothing
    # checked that.
    N = 64
    f = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
    # Unlike `f` in every slot, but still inside the PFC bounds: the point is to
    # leave a different pattern behind in the workspace, not to trip the check.
    junk = [1.9 - 1.8*(i-1)/N for i = 1:N]
    cases = vcat(uniform_schemes(fmin = 0.0, fmax = 2.0),
                 [("PFCNonUniform", PFCNonUniform(fill(0.05, N); fmin = 0.0, fmax = 2.0))])
    for (name, scheme) in cases
        c = name == "PFCNonUniform" ? 0.02 : 0.4
        fresh = similar(f)
        advect!(fresh, f, scheme, c, workspace(scheme, N))

        dirty = workspace(scheme, N)
        soiled = similar(f)
        advect!(soiled, junk, scheme, c, dirty)     # leave whatever it leaves
        advect!(soiled, f, scheme, c, dirty)
        @test soiled == fresh
    end
end

@testset "the kernels are type-stable" begin
    # The allocation gate notes a single boxed value appearing on some hosts and
    # not others. `@code_warntype` was clean when that was written; this keeps it
    # that way, and catches an inference regression before it shows up as
    # allocation on one CI runner only.
    N = 64
    src = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
    dst = similar(src)
    cases = vcat(uniform_schemes(fmin = 0.0, fmax = 2.0),
                 [("PFCNonUniform", PFCNonUniform(fill(0.05, N); fmin = 0.0, fmax = 2.0))])
    for (name, scheme) in cases
        c = name == "PFCNonUniform" ? 0.02 : 0.4
        @test (@inferred advect!(dst, src, scheme, c, workspace(scheme, N))) === dst
    end
end

@testset "element types" begin
    # Float32 data runs, and the result keeps its type. The workspaces do not:
    # `workspace(::SemiLagrangian, n)` hard-codes `Vector{Float64}`, so a Float32
    # line is prefiltered in double precision and narrowed on the way out.
    # Pinned as it stands rather than fixed -- making the workspaces generic is a
    # change to `src`, and worth doing on its own.
    N = 64
    g = Float32[1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
    for (name, scheme) in [("Upwind", Upwind()), ("LaxWendroff", LaxWendroff()),
                           ("Godunov", Godunov(PiecewiseLinear(), VanLeer())),
                           ("SemiLagrangian", SemiLagrangian(CubicSpline())),
                           ("PFC", PFC(fmin = 0.0f0, fmax = 2.0f0))]
        dst = similar(g)
        advect!(dst, g, scheme, 0.4f0, workspace(scheme, N))
        @test eltype(dst) === Float32
        @test all(isfinite, dst)
    end

    # PFCNonUniform does carry its element type through, workspace included.
    s32 = PFCNonUniform(fill(0.05f0, N); fmin = 0.0f0, fmax = 2.0f0)
    @test s32 isa PFCNonUniform{Float32}
    @test workspace(s32, N).accumulator isa Vector{Float32}

    # The semi-Lagrangian workspace is the one that does not.
    @test workspace(SemiLagrangian(CubicSpline()), N).buffer isa Vector{Float64}
end
