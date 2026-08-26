using Vasilek: LaxWendroff, Upwind, Godunov, SemiLagrangian, PFC

"""
    test_1d_advection_step(mod, Δx, Δt, v, f₀, f₁, exp_norm_dev, args...; kwargs...)

Advect `f₀` by one step of `mod` and compare against `f₁`.

Both `generate_solver` overloads are exercised — the one with the Courant number
baked in at construction, and the one taking it at call time. They must agree:
these two drifted apart in `PFC` once already.
"""
function test_1d_advection_step(mod, Δx, Δt, v, f₀, f₁, exp_norm_dev, args...; kwargs...)
    c = v*Δt/Δx
    name = nameof(mod)

    f_baked = similar(f₀)
    advect! = mod.generate_solver(f₀, f_baked, c, args...; kwargs...)
    advect!()
    dev_baked = Δx*norm(f_baked - f₁)
    println("$name $args $(values(values(kwargs))) (constant velocity): $dev_baked")
    @test dev_baked ≈ 0 atol=exp_norm_dev

    f_call = similar(f₀)
    advect! = mod.generate_solver(f₀, f_call, args...; kwargs...)
    advect!(c)
    dev_call = Δx*norm(f_call - f₁)
    println("$name $args $(values(values(kwargs))): $dev_call")
    @test dev_call ≈ 0 atol=exp_norm_dev

    # the two overloads must compute the same thing, bit for bit
    @test f_baked == f_call
end

@testset "Test 1D advection solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    v = 1.0

    f₀ = [1.0 + 0.01*sin(2π*i*Δx) for i = 0:99]
    f₁ = [1.0 + 0.01*sin(2π*(i*Δx - v*Δt)) for i = 0:99]

    test_1d_advection_step(Upwind, Δx, Δt,  v, f₀, f₁, 3e-7)
    test_1d_advection_step(Upwind, Δx, Δt, -v, f₁, f₀, 3e-7)

    test_1d_advection_step(LaxWendroff, Δx, Δt,  v, f₀, f₁, 1e-8)
    test_1d_advection_step(LaxWendroff, Δx, Δt, -v, f₁, f₀, 1e-8)

    for riemann_solver in (:Riemann_constant, :Riemann_linear)
        test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1e-6, riemann_solver)
        test_1d_advection_step(Godunov, Δx, Δt, -v, f₁, f₀, 1e-6, riemann_solver)
        if riemann_solver != :Riemann_constant
            for flux_limiter in (:VanLeer,)
                test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1e-6, riemann_solver; flux_limiter=flux_limiter)
                test_1d_advection_step(Godunov, Δx, Δt, -v, f₁, f₀, 1e-6, riemann_solver; flux_limiter=flux_limiter)
            end
        end
    end

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 3e-7; interpolation_order = :Linear)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₁, f₀, 3e-7; interpolation_order = :Linear)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 2e-9; interpolation_order = :Quadratic)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₁, f₀, 2e-9; interpolation_order = :Quadratic)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 2e-11; interpolation_order = :Cubic)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₁, f₀, 2e-11; interpolation_order = :Cubic)

    test_1d_advection_step(PFC, Δx, Δt,  v, f₀, f₁, 3e-8; fₘᵢₙ=0.0, fₘₐₓ=maximum(f₀))
    test_1d_advection_step(PFC, Δx, Δt, -v, f₁, f₀, 3e-8; fₘᵢₙ=0.0, fₘₐₓ=maximum(f₁))

    f₀ = [exp(-((i*Δx - 0.5)/0.15)^2) for i = 0:100]
    f₁ = [exp(-((i*Δx - v*Δt - 0.5)/0.15)^2) for i = 0:100]

    test_1d_advection_step(Upwind, Δx, Δt,  v, f₀, f₁, 1e-4)
    test_1d_advection_step(Upwind, Δx, Δt, -v, f₁, f₀, 1e-4)

    test_1d_advection_step(LaxWendroff, Δx, Δt,  v, f₀, f₁, 1e-5)
    test_1d_advection_step(LaxWendroff, Δx, Δt, -v, f₁, f₀, 1e-5)

    test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1e-4, :Riemann_constant)
    test_1d_advection_step(Godunov, Δx, Δt, -v, f₁, f₀, 1e-4, :Riemann_constant)

    test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1.1e-4, :Riemann_linear)
    test_1d_advection_step(Godunov, Δx, Δt, -v, f₁, f₀, 1.1e-4, :Riemann_linear)

    test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1.2e-4, :Riemann_linear; flux_limiter=:VanLeer)
    test_1d_advection_step(Godunov, Δx, Δt, -v, f₁, f₀, 1.2e-4, :Riemann_linear; flux_limiter=:VanLeer)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 3e-5; interpolation_order = :Linear)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₁, f₀, 3e-5; interpolation_order = :Linear)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 4e-7; interpolation_order = :Quadratic)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₁, f₀, 4e-7; interpolation_order = :Quadratic)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 4e-8; interpolation_order = :Cubic)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₁, f₀, 4e-8; interpolation_order = :Cubic)

    test_1d_advection_step(PFC, Δx, Δt,  v, f₀, f₁, 5e-6; fₘᵢₙ=0.0, fₘₐₓ=maximum(f₀))
    test_1d_advection_step(PFC, Δx, Δt, -v, f₁, f₀, 5e-6; fₘᵢₙ=0.0, fₘₐₓ=maximum(f₁))

    f₀ = zeros(Float64, 100)
    for i = 40:50
        f₀[i] = 1.0
    end
    f₁ = zeros(Float64, 100)
    for i = 41:50
        f₁[i] = 1.0
    end
    f₁[40] = 0.2
    f₁[51] = 0.8
    f₂ = zeros(Float64, 100)
    for i = 40:49
        f₂[i] = 1.0
    end
    f₂[39] = 0.8
    f₂[50] = 0.2

    test_1d_advection_step(Upwind, Δx, Δt,  v, f₀, f₁, 1e-18)
    test_1d_advection_step(Upwind, Δx, Δt, -v, f₀, f₂, 1e-18)

    test_1d_advection_step(LaxWendroff, Δx, Δt,  v, f₀, f₁, 1.6e-3)
    test_1d_advection_step(LaxWendroff, Δx, Δt, -v, f₀, f₂, 1.6e-3)

    test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1e-18, :Riemann_constant)
    test_1d_advection_step(Godunov, Δx, Δt, -v, f₀, f₂, 1e-18, :Riemann_constant)

    test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1e-2, :Riemann_linear)
    test_1d_advection_step(Godunov, Δx, Δt, -v, f₀, f₂, 1e-2, :Riemann_linear)

    test_1d_advection_step(Godunov, Δx, Δt,  v, f₀, f₁, 1e-2, :Riemann_linear; flux_limiter=:VanLeer)
    test_1d_advection_step(Godunov, Δx, Δt, -v, f₀, f₂, 1e-2, :Riemann_linear; flux_limiter=:VanLeer)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 5e-17; interpolation_order = :Linear)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₀, f₂, 5e-17; interpolation_order = :Linear)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 2e-3; interpolation_order = :Quadratic)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₀, f₂, 2e-3; interpolation_order = :Quadratic)

    test_1d_advection_step(SemiLagrangian, Δx, Δt,  v, f₀, f₁, 2e-3; interpolation_order = :Cubic)
    test_1d_advection_step(SemiLagrangian, Δx, Δt, -v, f₀, f₂, 2e-3; interpolation_order = :Cubic)

    test_1d_advection_step(PFC, Δx, Δt,  v, f₀, f₁, 1e-18; fₘᵢₙ=0.0, fₘₐₓ=maximum(f₀))
    test_1d_advection_step(PFC, Δx, Δt, -v, f₀, f₂, 1e-18; fₘᵢₙ=0.0, fₘₐₓ=maximum(f₀))
end
