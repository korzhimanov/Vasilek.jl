using Vasilek: Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC, PFCNonUniform,
               PoissonFourier1D, FDTD1D, BGK, Landau1P

"""
Allocation gate.

`@allocated` is deterministic on a given machine, unlike wall-clock timing on
a shared CI runner, so this can block a build without flapping. It catches the
regression class the refactors risk: a buffer that stops being reused, or a
closure capture that starts boxing.

**Why these are not asserted at an absolute zero.** On Julia 1.10 the optimizer
leaves a single boxed value (16 bytes) in one kernel or another depending on
the host it compiles for -- `SemiLagrangian` linear on all three CI runners,
`BGK` on the Windows runner, neither on the development machine, all at the
same patch version 1.10.12. `@code_warntype` is clean, so it is escape
analysis rather than an inference failure, and three attempts to remove it from
`SemiLagrangian` either did nothing or made it worse (64 bytes when the closure
was split, 32-40 kB when the interpolation object was hoisted). Julia 1.12
elides it everywhere.

So O(1) is asserted as **independence of N** plus a small ceiling. That still
fails on a reallocated buffer or a per-element temporary, which are the
regressions worth catching, and tolerates one constant box.

Everything runs inside a function. `@allocated` at global scope with
non-`const` globals measures the harness boxing its own arguments -- while
preparing this file a stray `c*Δx` built from two globals reported 16 bytes
for a solver that allocates nothing.
"""
const ALLOC_N = 1000
const ALLOC_CEILING = 1024

alloc_data(n) = [1.0 + 0.5*sin(2π*i/n) for i = 0:n-1]

"Bytes allocated by one warmed-up step of the solver built by `mk(src, dst)`."
function step_bytes(mk, n, c)
    src = alloc_data(n)
    dst = similar(src)
    step! = mk(src, dst)
    step!(c)
    step!(c)
    return @allocated step!(c)
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
            ("Upwind",           (a,b) -> Upwind.generate_solver(a,b)),
            ("LaxWendroff",      (a,b) -> LaxWendroff.generate_solver(a,b)),
            ("Godunov constant", (a,b) -> Godunov.generate_solver(a,b,:Riemann_constant)),
            ("Godunov linear",   (a,b) -> Godunov.generate_solver(a,b,:Riemann_linear)),
            ("Godunov VanLeer",  (a,b) -> Godunov.generate_solver(a,b,:Riemann_linear; flux_limiter=:VanLeer)),
            ("PFC",              (a,b) -> PFC.generate_solver(a,b;fₘᵢₙ=0.0,fₘₐₓ=2.0)),
            ("SemiLagrangian linear",
                                 (a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=:Linear)),
        ]
        for (name, mk) in cases
            assert_constant(name, n -> step_bytes(mk, n, 0.4))
        end
    end

    @testset "SemiLagrangian spline prefilter" begin
        # These genuinely scale with N: the allocation is inside Interpolations'
        # periodic prefilter, not in this package. copyto! into the preallocated
        # buffer is free, and `extrapolate` and the evaluation loop add nothing;
        # `interpolate!` accounts for all of it. Hoisting the buffer took
        # quadratic and cubic from 172_868 and 172_884 to about 100_700 at
        # N = 1000, so roughly 100 bytes per element. Bounded so a regression is
        # still visible.
        for order in (:Quadratic, :Cubic)
            bytes = step_bytes((a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=order),
                               ALLOC_N, 0.4)
            println("  ", rpad("SemiLagrangian $order", 24), bytes)
            @test bytes < 110*ALLOC_N
        end
    end

    @testset "non-uniform advection, fields and collisions are O(1)" begin
        function nonuniform_bytes(n)
            Δ = fill(0.01, n)
            f = alloc_data(n)
            step! = PFCNonUniform.make_advect_1D!(Δ; fₘᵢₙ = 0.0, fₘₐₓ = 2.0)
            step!(f, 0.004); step!(f, 0.004)
            return @allocated step!(f, 0.004)
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

        function bgk_bytes(n)
            v = collect(range(-4, 4; length = n))
            f = @. exp(-v^2)
            step! = BGK.generate_solver(f, v, 0.1, 1e-2)
            step!(); step!()
            return @allocated step!()
        end

        function landau_bytes(n)
            v = collect(range(-4, 4; length = n))
            src = @. exp(-v^2)
            f = similar(src)
            step! = Landau1P.generate_solver(src, f, v, 0.1, 1e-2)
            step!(); step!()
            return @allocated step!()
        end

        assert_constant("PFCNonUniform", nonuniform_bytes)
        assert_constant("PoissonFourier1D", poisson_bytes)
        assert_constant("FDTD1D", fdtd_bytes)
        assert_constant("BGK", bgk_bytes)
        # Landau1P is O(N^2) in time; a smaller pair keeps the test quick.
        assert_constant("Landau1P", n -> landau_bytes(n ÷ 10))
    end
end
