using Vasilek
using Vasilek.Collisions: Landau1P, collide!, workspace as collision_workspace
using Vasilek: PoissonFourier1D, FDTD1D

include(joinpath(@__DIR__, "scheme_cases.jl"))

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
    println("  ", rpad(name, 24), "N=", ALLOC_N, ": ", small, "   N=", 4*ALLOC_N, ": ", large)
    @test small == large
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
            println("  ", rpad("SemiLagrangian $(typeof(spline).name.name)", 24), bytes)
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
            solve! = PoissonFourier1D.generate_solver(ρ, 0.01)
            solve!(e, ρ); solve!(e, ρ)
            return @allocated solve!(e, ρ)
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
            ws = collision_workspace(op, n)
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
end
