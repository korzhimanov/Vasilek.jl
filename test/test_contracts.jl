using Vasilek

@isdefined(march!) || include(joinpath(@__DIR__, "scheme_cases.jl"))

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
            @test_throws DomainError march!(over, PFC(fmin = 0.0, fmax = 2.0), c, 1)
            @test_throws DomainError march!(under, PFC(fmin = 0.0, fmax = 2.0), c, 1)
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
        # and the flag is a Bool: `PFC{Float64, 3}` used to construct, and
        # then stop its first step with a TypeError from `if Checked`
        @test_throws ArgumentError PFC{Float64, 3}(0.0, 2.0)
        @test PFC{Float64, true}(0.0, 2.0) === PFC(fmin = 0.0, fmax = 2.0)
    end

    @testset "a constant line at fmax stays constant, to the bit" begin
        # Every face of a constant line carries the same flux, c*f, and each
        # cell takes the difference of its two, which is zero. It used to add
        # the inflow and subtract the outflow in turn, and the two roundings
        # left a cell an ulp above fmax, which the check above then refuses:
        # at 42 of the 401 Courant numbers below from 0.7, and at 108 from
        # 1/√(2π).
        for v in (0.7, 1/sqrt(2π))
            line = fill(v, N)
            scheme = PFC(fmin = 0.0, fmax = v)
            ws = workspace(scheme, N)
            @test all(c -> advect!(similar(line), line, scheme, c, ws) == line,
                      range(-1, 1; length = 401))
        end
    end
end

# `PFCNonUniform` takes the same bounds and builds them into the same limiter,
# and checked none of them: the verification harness bounded every run at
# `fmax = 1`, an equilibrium peaking at 1.79 went through it, and it ended 43.6%
# of its peak away from itself by t = 50 without an error. The same tests as
# `PFC`'s above, on the same data and at the same Courant number -- the fourth
# argument here is a displacement, and 0.02 over cells of 0.05 is 0.4.
@testset "PFCNonUniform bounds check" begin
    N = 64
    f = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]      # ranges over [0.5, 1.5]
    Δx = fill(0.05, N)
    scheme(; kw...) = PFCNonUniform(Δx; fmin = 0.0, fmax = 2.0, kw...)

    @testset "the check fires, in both branches" begin
        over = copy(f);  over[10] = 5.0
        under = copy(f); under[10] = -1.0
        for α in (0.02, -0.02)
            @test_throws DomainError march!(over, scheme(), α, 1)
            @test_throws DomainError march!(under, scheme(), α, 1)
        end

        # Touching the bounds exactly is legal: the assertions are `≤`.
        exact = copy(f); exact[1] = 0.0; exact[2] = 2.0
        @test march!(exact, scheme(), 0.02, 1) isa Vector{Float64}
    end

    @testset "and says what PFC says" begin
        # One function checks both schemes, so the same data draws the same
        # message from either: the bound, and the extremum that broke it.
        over = copy(f);  over[10] = 5.0
        under = copy(f); under[10] = -1.0
        message(s, c, g) = try march!(g, s, c, 1); "" catch e; e.msg end
        @test message(scheme(), 0.02, over) == message(PFC(fmin = 0.0, fmax = 2.0), 0.4, over) ==
              "fmax = 2.0 is below maximum(src) = 5.0: the PFC limiter is built on [fmin, fmax]"
        @test message(scheme(), 0.02, under) == message(PFC(fmin = 0.0, fmax = 2.0), 0.4, under) ==
              "fmin = 0.0 exceeds minimum(src) = -1.0: the PFC limiter is built on [fmin, fmax]"
    end

    @testset "checked = false compiles it away without changing the answer" begin
        for α in (0.02, -0.02)
            @test march!(f, scheme(checked = true), α, 4) == march!(f, scheme(checked = false), α, 4)
        end

        # On data outside the bounds it does not throw, and the answer it gives
        # instead does not look wrong: one step from a peak of 5 against
        # `fmax = 2` lands 0.048 from the same step bounded at 5 -- 0.37 after
        # four -- and conserves mass to round-off, so a mass check downstream
        # would pass it.
        over = copy(f); over[10] = 5.0
        wrong = march!(over, scheme(checked = false), 0.02, 1)
        right = march!(over, PFCNonUniform(Δx; fmin = 0.0, fmax = 5.0), 0.02, 1)
        @test maximum(abs, wrong .- right) > 0.01
        @test sum(wrong) ≈ sum(over) rtol = 1e-13
    end

    @testset "the flag is a type parameter, and construction still infers" begin
        @test scheme()                  isa PFCNonUniform{Float64, true}
        @test scheme(checked = false)   isa PFCNonUniform{Float64, false}
        # The one-parameter spelling in existing code matches either.
        @test scheme()                  isa PFCNonUniform{Float64}
        @test scheme(checked = false)   isa PFCNonUniform{Float64}
        # The element type is the grid's, as before, and the bounds follow it.
        @test PFCNonUniform(fill(0.05f0, N); fmin = 0, fmax = 2) isa PFCNonUniform{Float32, true}
        # `checked` is a value, and it has to reach the type. Inference did not
        # carry the default through the constructor on its own, so the call
        # below inferred as `PFCNonUniform{Float64}` -- concrete before the flag
        # existed, abstract after -- until `@constprop :aggressive`.
        @test (@inferred PFCNonUniform(Δx; fmin = 0.0, fmax = 2.0)) isa PFCNonUniform{Float64, true}
        # and the flag is a Bool here too
        @test_throws ArgumentError PFCNonUniform{Float64, 3}(Δx, fill(2.0, N), 0.05, 0.0, 2.0)
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

    @testset "so is a view sharing src's memory" begin
        # `===` alone let these through: `view(f, :)` as dest was 0.047 off for
        # Upwind with no error. `Base.mightalias` sees the shared memory.
        for (name, scheme) in all_schemes
            a = copy(f)
            c = name == "PFCNonUniform" ? 0.02 : 0.4
            @test_throws ArgumentError advect!(view(a, :), a, scheme, c, workspace(scheme, N))
            @test_throws ArgumentError advect!(a, view(a, :), scheme, c, workspace(scheme, N))
        end
        # disjoint halves of one array are not aliases
        buf = repeat(f, 2)
        @test advect!(view(buf, 1:N), view(buf, N+1:2N), Upwind(), 0.4) == advect!(similar(f), f, Upwind(), 0.4)
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

    @testset "so is a workspace for a scheme that takes none" begin
        # Upwind and the rest took any workspace at all and ignored it, so a
        # sweep that handed one scheme another's workspace went unnoticed until
        # a scheme that reads its workspace was handed the wrong one.
        foreign = (workspace(SemiLagrangian(CubicSpline()), N),
                   workspace(PFCNonUniform(fill(0.05, N); fmin = 0.0, fmax = 2.0), N))
        for (name, scheme) in uniform_schemes(fmin = 0.0, fmax = 2.0), ws in foreign
            startswith(name, "SemiLagrangian") && continue
            @test_throws ArgumentError advect!(similar(f), f, scheme, 0.4, ws)
        end
    end

    @testset "and a workspace sharing memory with dest or src" begin
        # Measured before the check: a dest sharing a spline buffer came back
        # 0.013 (linear) and 0.014 (cubic) off; a src sharing the cubic one was
        # overwritten with its coefficients, 8.0e-4 from what was passed; and a
        # src sharing PFCNonUniform's accumulator came back 0.21 off, 0.81 on a
        # 1:2 grid. A dest sharing the accumulator happened to come out right,
        # and is refused with the rest: the workspace is the scheme's own.
        for spline in (LinearSpline(), CubicSpline())
            s = SemiLagrangian(spline)
            ws = workspace(s, N)
            buf = view(ws.buffer, 1:N)
            copyto!(buf, f)
            @test_throws ArgumentError advect!(buf, f, s, 0.4, ws)
            @test_throws ArgumentError advect!(similar(f), buf, s, 0.4, ws)
        end
        s = PFCNonUniform(fill(0.05, N); fmin = 0.0, fmax = 2.0)
        ws = workspace(s, N)
        copyto!(ws.accumulator, f)
        @test_throws ArgumentError advect!(ws.accumulator, f, s, 0.02, ws)
        @test_throws ArgumentError advect!(similar(f), ws.accumulator, s, 0.02, ws)
        # A workspace type with no method of its own is an error, where it was
        # taken as sharing nothing: a new scheme's workspace cannot skip the check.
        @test_throws MethodError Vasilek.Advection._scratch_aliases((buffer = f,), f, similar(f))
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

# The fourth check `advect!` runs: the Courant bound. Every scheme but
# SemiLagrangian is explicit and unstable past |c| = 1, and one step there
# looks right -- LaxWendroff at c = 1.2 is 5.2e-6 from the exact shift of a
# smooth profile, against 1.7e-6 at 0.9 -- where a hundred steps are 9.1e10
# from it. The measurements are in the `_validate_courant` docstring.
@testset "advect! refuses a Courant number the scheme cannot take" begin
    N = 64
    f = [1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
    bounded = [(name, s) for (name, s) in uniform_schemes(fmin = 0.0, fmax = 2.0)
               if !startswith(name, "SemiLagrangian")]

    @testset "past |c| = 1 is refused, both ways, before anything is written" begin
        # The first float above one, a plain case, and the free-streaming row
        # that found this -- v = 6, Δt = 0.02, Δx = 4π/128 -- which PFC took
        # silently until its own bounds assertion fired 218 steps later.
        free_streaming_row = 6*0.02/(4π/128)
        for (name, scheme) in bounded, c in (nextfloat(1.0), free_streaming_row, 2.0),
                sgn in (1, -1)
            dst = fill(-7.0, N)
            @test_throws DomainError advect!(dst, f, scheme, sgn*c, workspace(scheme, N))
            @test all(==(-7.0), dst)
        end
    end

    @testset "c = ±1 exactly is accepted, and is a one-cell shift" begin
        # The bound is `≤`. Measured at c = ±1: Upwind exact; LaxWendroff, the
        # three Godunov and PFC at most one ulp (2.2e-16) from circshift.
        # Godunov(PiecewiseLinear()) is a shift here because its slope term
        # carries (1 − |c|), which vanishes. Before it did, the scheme was 2.4e-3
        # off without a limiter and 2.7e-3 with VanLeer, and this test excused it.
        for (name, scheme) in bounded, c in (1.0, -1.0)
            out = march!(f, scheme, c, 1)
            @test maximum(abs, out .- circshift(f, Int(c))) ≤ 2eps()
        end
    end

    @testset "NaN is refused" begin
        # No comparison with NaN holds, so the check refuses it without a case
        # of its own. It used to go through and come back NaN everywhere.
        for (name, scheme) in bounded
            @test_throws DomainError advect!(similar(f), f, scheme, NaN, workspace(scheme, N))
        end
    end

    @testset "checked = false does not switch it off" begin
        # `checked` compiles away PFC's minimum/maximum pass. The Courant check
        # is one comparison and stays: a PFC run past the limit is exactly what
        # that pass used to be the only thing catching.
        @test_throws DomainError advect!(similar(f), f,
                                         PFC(fmin = 0.0, fmax = 2.0, checked = false), 1.2)
    end

    @testset "SemiLagrangian refuses only a non-finite c" begin
        # It used to take NaN and Inf too, and return NaN everywhere.
        for spline in (LinearSpline(), CubicSpline()), c in (NaN, Inf, -Inf)
            s = SemiLagrangian(spline)
            @test_throws DomainError advect!(similar(f), f, s, c, workspace(s, N))
        end
    end

    @testset "SemiLagrangian has no Courant limit, and is not checked" begin
        # Its accuracy past |c| = 1 is test_symmetry's business (c = 3.7, and
        # c ± N); here, only that the check does not reach it.
        for spline in (LinearSpline(), QuadraticSpline(), CubicSpline()), c in (1.2, -2.0, 3.7)
            s = SemiLagrangian(spline)
            @test advect!(similar(f), f, s, c, workspace(s, N)) isa Vector{Float64}
        end
    end

    @testset "PFCNonUniform is held to its narrowest cell" begin
        # Its fourth argument is a displacement, and every cell gives up its
        # flux alone, so the bound is the narrowest width however wide the rest
        # are. On a 1:2 grid a step of 1.1 narrow widths -- 0.55 of the wide one
        # -- took data with an exact zero to -7.6e-6.
        s = PFCNonUniform(vcat(fill(0.1, 16), fill(0.05, 32), fill(0.1, 16));
                          fmin = 0.0, fmax = 2.0)
        ws = workspace(s, N)
        for α in (0.05, -0.05)
            @test advect!(similar(f), f, s, α, ws) isa Vector{Float64}
        end
        for α in (nextfloat(0.05), -nextfloat(0.05), 0.055, 0.1)
            @test_throws DomainError advect!(similar(f), f, s, α, ws)
        end
    end

    @testset "the error carries the value and names the scheme" begin
        err = try advect!(similar(f), f, LaxWendroff(), 1.2) catch e; e end
        @test err isa DomainError && err.val == 1.2
        @test occursin("LaxWendroff", sprint(showerror, err))
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
    # Float32 data runs and keeps its type, and so do the workspaces, built
    # with `workspace(scheme, n, Float32)` or by the four-argument `advect!`
    # from `src`. `PFCNonUniform` used to throw a MethodError on Float32 data
    # (its limiter was annotated `::Float64`), and the spline workspace was
    # Float64 whatever the data.
    N = 64
    g = Float32[1.0 + 0.5*sin(2π*(i-1)/N) for i = 1:N]
    for (name, scheme, c) in [("Upwind", Upwind(), 0.4f0), ("LaxWendroff", LaxWendroff(), 0.4f0),
                              ("Godunov", Godunov(PiecewiseLinear(), VanLeer()), 0.4f0),
                              ("Godunov Superbee", Godunov(PiecewiseLinear(), Superbee()), 0.4f0),
                              ("SemiLagrangian", SemiLagrangian(CubicSpline()), 0.4f0),
                              ("SemiLagrangian linear", SemiLagrangian(LinearSpline()), 0.4f0),
                              ("PFC", PFC(fmin = 0.0f0, fmax = 2.0f0), 0.4f0),
                              ("PFCNonUniform", PFCNonUniform(fill(0.05f0, N); fmin = 0.0f0, fmax = 2.0f0), 0.02f0)]
        dst = similar(g)
        @test (@inferred advect!(dst, g, scheme, c, workspace(scheme, N, Float32))) === dst
        @test all(isfinite, dst)
        # and the answer is the Float64 one to single precision
        ref = advect!(similar(g, Float64), Float64.(g), scheme, Float64(c))
        @test maximum(abs, dst .- ref) < 1e-5
        @test eltype(advect!(similar(g), g, scheme, c)) === Float32
    end

    @test workspace(SemiLagrangian(CubicSpline()), N, Float32).buffer isa Vector{Float32}
    @test workspace(SemiLagrangian(CubicSpline()), N).buffer isa Vector{Float64}
    s32 = PFCNonUniform(fill(0.05f0, N); fmin = 0.0f0, fmax = 2.0f0)
    @test workspace(s32, N).accumulator isa Vector{Float32}

    # Bounds narrower than they were given round outward, so a Float32 grid
    # keeps the Float64 data they admit: `fmax = maximum(f)` rounded below the
    # peak as often as not, and the first checked step refused the data's own.
    peak = 0.4029317030250849
    @test Float32(peak) < peak
    s = PFCNonUniform(fill(0.05f0, N); fmin = -peak, fmax = peak)
    @test s isa PFCNonUniform{Float32, true}
    @test s.fmin ≤ -peak && s.fmax ≥ peak
    line = [peak*sin(2π*(i-1)/N) for i in 1:N]
    line[N÷4 + 1] = peak
    @test all(isfinite, advect!(similar(line), line, s, 0.02f0))

    # The driver's cell widths keep a Float32 grid's type, where a `0.5` in
    # them used to make Float64 of it.
    widths = Vasilek.VlasovPoisson1D1V.cell_widths(Float32[0, 1, 3, 4])
    @test eltype(widths) === Float32
    @test widths == Float32[1, 1.5, 1.5, 1]
    # and they are the widths `BGK` weighs its moments by, a second copy of the
    # formula; the two had drifted once, `0.5*` in one and `/2` in the other
    σ = range(-1, 1; length = 41)
    for z in (@.(6sinh(2σ)/sinh(2)), Float32.(@.(6sinh(2σ)/sinh(2))), cumsum(1 .+ sin.(1:20).^2))
        w = Vasilek.VlasovPoisson1D1V.cell_widths(z)
        @test w == [Vasilek.Collisions._width(z, i) for i in eachindex(z)]
        @test eltype(w) === eltype(z)
    end
end

@testset "constructors refuse what the schemes cannot use" begin
    @test_throws ArgumentError PFC(fmin = 2.0, fmax = 1.0)
    @test_throws ArgumentError PFCNonUniform(fill(0.1, 8); fmin = 2.0, fmax = 1.0)
    @test_throws ArgumentError PFCNonUniform([0.1, 0.0, 0.1, 0.1]; fmin = 0.0, fmax = 1.0)
    @test_throws ArgumentError PFCNonUniform([0.1, -0.1, 0.1, 0.1]; fmin = 0.0, fmax = 1.0)
    @test_throws ArgumentError PFCNonUniform([0.1, NaN, 0.1, 0.1]; fmin = 0.0, fmax = 1.0)
    @test_throws ArgumentError PFCNonUniform([0.1, 0.1]; fmin = 0.0, fmax = 1.0)
    # and so do the positional forms every scheme is built through, which
    # checked the flag alone: `PFC{Float64, false}(2.0, 0.0)` built, and a
    # PFCNonUniform told its narrowest cell was 1.0, over cells of 0.05, took a
    # 20-cell step from data in [0.5, 1.5] to [-1.02, 3.02] without an error
    Δx = fill(0.05, 40)
    @test_throws ArgumentError PFC{Float64, false}(2.0, 0.0)
    @test_throws ArgumentError PFCNonUniform{Float64, true}(Δx, fill(2.0, 40), 1.0, 0.0, 2.0)
    @test_throws ArgumentError PFCNonUniform{Float64, true}(Δx, fill(2.0, 40), 0.05, 2.0, 0.0)
    @test_throws ArgumentError PFCNonUniform{Float64, true}([0.1, 0.1], [2.0, 2.0], 0.1, 0.0, 1.0)
    @test_throws ArgumentError PFCNonUniform{Float64, true}(-Δx, fill(2.0, 40), -0.05, 0.0, 1.0)
    @test_throws DimensionMismatch PFCNonUniform{Float64, true}(Δx, fill(2.0, 39), 0.05, 0.0, 2.0)
    @test PFCNonUniform{Float64, true}(Δx, fill(2.0, 40), 0.05, 0.0, 2.0).Δxmin == 0.05
    # a piecewise-constant reconstruction has no slope to limit
    @test_throws ArgumentError Godunov(PiecewiseConstant(), VanLeer())
    @test_throws ArgumentError Godunov(PiecewiseConstant(), Superbee())
    @test Godunov(PiecewiseConstant()) isa Godunov{PiecewiseConstant, NoLimiter}
    @test Godunov(PiecewiseConstant(), NoLimiter()) isa Godunov{PiecewiseConstant, NoLimiter}
end
