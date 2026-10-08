using Vasilek
using Vasilek.Collisions: Landau1P, collide!
using Vasilek: PoissonFourier1D, FDTD1D

"""
Allocation gate.

`@allocated` is deterministic on a given machine, unlike wall-clock timing on a
shared CI runner, so this can block a build without flapping. It catches the
regression class the refactors risk: a buffer that stops being reused, or a
capture that starts boxing.

**Why these are not asserted at an absolute zero.** On Julia 1.10 the optimizer
leaves a single boxed value (16 bytes) in one kernel or another depending on
the host it compiles for -- `SemiLagrangian` linear on all three CI runners,
`BGK` on the Windows runner, neither on the development machine, all at the
same patch version 1.10.12. `@code_warntype` is clean, so it is escape analysis
rather than an inference failure, and Julia 1.12 elides it everywhere.

So O(1) is asserted as **independence of N** plus a small ceiling. That still
fails on a reallocated buffer or a per-element temporary, which are the
regressions worth catching, and tolerates one constant box.

Everything runs inside a function. `@allocated` at global scope with non-`const`
globals measures the harness boxing its own arguments.
"""
const ALLOC_N = 1000
const ALLOC_CEILING = 1024

alloc_data(n) = [1.0 + 0.5*sin(2π*i/n) for i = 0:n-1]

"Bytes allocated by one warmed-up step of `scheme` at size `n`."
function step_bytes(scheme, n, c)
    src = alloc_data(n)
    dst = similar(src)
    ws = workspace(scheme, n)
    advect!(dst, src, scheme, c, ws)
    advect!(dst, src, scheme, c, ws)
    return @allocated advect!(dst, src, scheme, c, ws)
end

"Assert the allocation count is O(1) in the problem size."
function assert_constant(name, bytes_at)
    small = bytes_at(ALLOC_N)
    large = bytes_at(4*ALLOC_N)
    println("  ", rpad(name, 32), "N=", ALLOC_N, ": ", small, "   N=", 4*ALLOC_N, ": ", large)
    # Independent of N to within 64 bytes, four of the 16-byte boxes described
    # above: an exact equality would also fail on escape analysis deciding
    # differently at the two sizes, which is not the regression this guards
    # against, while a reallocated buffer or a per-element temporary at N = 1000
    # is kilobytes.
    @test abs(large - small) ≤ 64
    @test small < ALLOC_CEILING
end

@testset "Allocations" begin
    @testset "advection kernels are O(1)" begin
        cases = [
            ("Upwind",                Upwind()),
            ("LaxWendroff",           LaxWendroff()),
            ("Godunov constant",      Godunov(PiecewiseConstant())),
            ("Godunov linear",        Godunov(PiecewiseLinear())),
            ("Godunov VanLeer",       Godunov(PiecewiseLinear(), VanLeer())),
            ("Godunov Superbee",      Godunov(PiecewiseLinear(), Superbee())),
            ("PFC",                   PFC(fmin = 0.0, fmax = 2.0)),
            ("SemiLagrangian linear", SemiLagrangian(LinearSpline())),
        ]
        for (name, scheme) in cases
            assert_constant(name, n -> step_bytes(scheme, n, 0.4))
        end
    end

    @testset "SemiLagrangian spline prefilter" begin
        # These genuinely scale with N: the allocation is inside the periodic
        # prefilter in Interpolations, not in this package. copyto! into the
        # preallocated buffer is free, and `extrapolate` and the evaluation loop
        # add nothing; `interpolate!` accounts for all of it. About 100 bytes
        # per element. Bounded so a regression is still visible.
        for spline in (QuadraticSpline(), CubicSpline())
            bytes = step_bytes(SemiLagrangian(spline), ALLOC_N, 0.4)
            println("  ", rpad("SemiLagrangian $(typeof(spline).name.name)", 32), bytes)
            @test bytes < 110*ALLOC_N
        end
    end

    @testset "non-uniform advection, fields and collisions are O(1)" begin
        function nonuniform_bytes(n)
            scheme = PFCNonUniform(fill(0.01, n); fmin = 0.0, fmax = 2.0)
            src = alloc_data(n)
            dst = similar(src)
            ws = workspace(scheme, n)
            advect!(dst, src, scheme, 0.004, ws); advect!(dst, src, scheme, 0.004, ws)
            return @allocated advect!(dst, src, scheme, 0.004, ws)
        end

        function poisson_bytes(n)
            ρ = [sin(2π*i/n) for i = 0:n-1]
            e = similar(ρ)
            p = PoissonFourier1D.PoissonFFT1D(n, 0.01); ws = workspace(p)
            PoissonFourier1D.solve!(e, ρ, p, ws); PoissonFourier1D.solve!(e, ρ, p, ws)
            return @allocated PoissonFourier1D.solve!(e, ρ, p, ws)
        end

        function fdtd_bytes(n)
            Δx = 0.01; Δt = 0.8*Δx
            mesh = FDTD1D.YeeMesh1D{Float64}(n)
            pulse = (y = (t,x) -> 0.0, z = (t,x) -> 0.0)
            advance! = FDTD1D.make_advance_fields(mesh, Δt/Δx, pulse, Δt, Δx, 0.0,
                                                  FDTD1D.PML(; N = 0, σ_max = 1.0, Δx = Δx, Δt = Δt))
            j = (y = zeros(n + 1), z = zeros(n + 1))
            advance!(0.0, j); advance!(0.0, j)
            return @allocated advance!(0.0, j)
        end

        function collision_bytes(op, n)
            v = collect(range(-4, 4; length = n))
            src = @. exp(-v^2)
            dst = similar(src)
            ws = workspace(op, n)
            collide!(dst, src, op, v, 0.1, ws); collide!(dst, src, op, v, 0.1, ws)
            return @allocated collide!(dst, src, op, v, 0.1, ws)
        end

        assert_constant("PFCNonUniform", nonuniform_bytes)
        assert_constant("PoissonFourier1D", poisson_bytes)
        assert_constant("FDTD1D", fdtd_bytes)
        assert_constant("BGK", n -> collision_bytes(BGK(1e-2), n))
        # Landau1P is O(N^2) in time; a smaller pair keeps the test quick.
        assert_constant("Landau1P", n -> collision_bytes(Landau1P(1e-2), n ÷ 10))
    end

    @testset "the driver's step is O(1)" begin
        # `vlasov_poisson` allocates its buffers once, before the first step. A
        # step is the difference between a 12-step and a 2-step run, over 10,
        # which cancels the setup; what remains is the histories' growth, eight
        # entries and one mode, 80 bytes a step in Float64, and anything the step
        # itself allocates. Measured 80 on Julia 1.10 and 1.13 at both sizes:
        # the step allocates nothing. It allocated 84 KiB a step at 64 × 161,
        # more than `f`'s 80.5 KiB, for the density summed from a temporary.
        function driver_bytes(Nx, Nv, steps; kw...)
            k = 0.5
            x = collect(range(2π/k/Nx; step = 2π/k/Nx, length = Nx))
            v = collect(range(-6.0, 6.0; length = Nv))
            f₀ = [exp(-u^2/2)/sqrt(2π)*(1 + 0.01cos(k*y)) for u in v, y in x]
            t = collect(0.0:0.05:0.05*steps)
            vlasov_poisson(x, v, f₀, t; modes = (k,), invariants = true, kw...)
            return @allocated vlasov_poisson(x, v, f₀, t; modes = (k,), invariants = true, kw...)
        end
        per_step(Nx, Nv; kw...) = (driver_bytes(Nx, Nv, 12; kw...) - driver_bytes(Nx, Nv, 2; kw...))/10

        small, large = per_step(64, 81), per_step(128, 161)
        println("  ", rpad("vlasov_poisson step", 32), "64 × 81: ", small, "   128 × 161: ", large)
        # Independent of the grid to four 16-byte boxes, and at most 112 bytes
        # above the histories' 80: a temporary the size of a line is 648 bytes
        # here, one the size of `f` 40.5 KiB.
        @test abs(large - small) ≤ 64
        @test small ≤ 192

        # With `BGK` the step is not held to that: on the Windows runner under
        # Julia 1.10 `collide!` leaves its 16-byte box (see the top of this
        # file), twice a line. It is held to a few boxes a line, at a fixed
        # number of lines, whatever the length of each.
        c81, c161 = per_step(64, 81; collisions = BGK(1.0)), per_step(64, 161; collisions = BGK(1.0))
        println("  ", rpad("vlasov_poisson step, BGK", 32), "64 × 81: ", c81, "   64 × 161: ", c161)
        @test abs(c161 - c81) ≤ 64
        @test c81 ≤ 192 + 64*64
    end
end
