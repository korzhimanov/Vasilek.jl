using Vasilek: Upwind, LaxWendroff, Godunov, SemiLagrangian, PFC, PFCNonUniform,
               PoissonFourier1D, FDTD1D, BGK, Landau1P

"""
Allocation gate.

`@allocated` is deterministic, unlike wall-clock timing on a shared CI runner,
so this can block a build without flapping. It catches exactly the class of
regression the upcoming refactors risk: a buffer that stops being reused, or a
closure capture that starts boxing.

Everything runs inside a function. `@allocated` at global scope with
non-`const` globals measures the harness boxing its own arguments — while
preparing this file a stray `c*Δx` built from two globals reported 16 bytes
for a solver that actually allocates nothing.
"""
const ALLOC_N = 1000

alloc_data() = [1.0 + 0.5*sin(2π*i*0.01) for i = 0:ALLOC_N-1]

function step_allocations(mk, c)
    src = alloc_data()
    dst = similar(src)
    step! = mk(src, dst)
    step!(c)                       # compile and warm every buffer
    step!(c)
    return @allocated step!(c)
end

@testset "Allocations" begin
    @testset "advection kernels allocate nothing" begin
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
            bytes = step_allocations(mk, 0.4)
            println("  ", rpad(name, 24), bytes)
            @test bytes == 0
        end
    end

    @testset "SemiLagrangian spline prefilter" begin
        # The remaining allocation is inside Interpolations' periodic B-spline
        # prefilter, not in this package: copyto! into the preallocated buffer
        # is free, `extrapolate` and the evaluation loop add nothing, and
        # `interpolate!` alone accounts for all of it. Hoisting the buffer took
        # quadratic and cubic from 172_868 and 172_884 down to about 100_700.
        #
        # Bounded rather than zero, so a regression is still visible.
        for order in (:Quadratic, :Cubic)
            bytes = step_allocations((a,b) -> SemiLagrangian.generate_solver(a,b;interpolation_order=order), 0.4)
            println("  ", rpad("SemiLagrangian $order", 24), bytes)
            @test bytes < 110_000
        end
    end

    @testset "non-uniform advection, fields and collisions" begin
        function nonuniform_bytes()
            Δ = fill(0.01, ALLOC_N)
            f = alloc_data()
            step! = PFCNonUniform.make_advect_1D!(Δ; fₘᵢₙ = 0.0, fₘₐₓ = 2.0)
            step!(f, 0.004); step!(f, 0.004)
            return @allocated step!(f, 0.004)
        end

        function poisson_bytes()
            ρ = [sin(2π*i*0.01) for i = 0:ALLOC_N-1]
            e = similar(ρ)
            solve! = PoissonFourier1D.generate_solver(ρ, 0.01)
            solve!(e, ρ); solve!(e, ρ)
            return @allocated solve!(e, ρ)
        end

        function fdtd_bytes()
            Δx = 0.01; Δt = 0.8*Δx
            mesh = FDTD1D.YeeMesh1D{Float64}(ALLOC_N)
            pulse = (y = (t,x) -> 0.0, z = (t,x) -> 0.0)
            advance! = FDTD1D.make_advance_fields(mesh, Δt/Δx, pulse, Δt, Δx, 0.0,
                                                  FDTD1D.PML(0, 1.0, Δx, Δt))
            j = (y = zeros(ALLOC_N + 1), z = zeros(ALLOC_N + 1))
            advance!(0.0, j); advance!(0.0, j)
            return @allocated advance!(0.0, j)
        end

        function bgk_bytes()
            v = collect(-4:0.01:4)
            f = @. exp(-v^2)
            step! = BGK.generate_solver(f, v, 0.1, 1e-2)
            step!(); step!()
            return @allocated step!()
        end

        function landau_bytes()
            v = collect(-4:0.1:4)
            src = @. exp(-v^2)
            f = similar(src)
            step! = Landau1P.generate_solver(src, f, v, 0.1, 1e-2)
            step!(); step!()
            return @allocated step!()
        end

        for (name, fn) in (("PFCNonUniform", nonuniform_bytes),
                           ("PoissonFourier1D", poisson_bytes),
                           ("FDTD1D", fdtd_bytes),
                           ("BGK", bgk_bytes),
                           ("Landau1P", landau_bytes))
            bytes = fn()
            println("  ", rpad(name, 24), bytes)
            @test bytes == 0
        end
    end
end
