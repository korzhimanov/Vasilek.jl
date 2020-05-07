const MODULES = (:LaxWendroff,
                 :Upwind,
                 :Godunov)

for mod_name in MODULES
    include(joinpath(dirname(@__FILE__),"..","..","src","VlasovSolver","$mod_name.jl"))
end

args_var = Dict(
    :Godunov => (:Riemann_constant, :Riemann_linear)
)

function test_1d_advection_step(mod_name, Δx, Δt, v, f₀, f₁, exp_norm_dev, args...; plot_needed=false)
    if plot_needed
        plt = plot(f₀, label = "initial")
        plot!(f₁, label = "expected")
    end

    f = similar(f₀)

    advect! = eval(:($mod_name)).generate_solver(f₀, f, v*Δt/Δx, args...)
    advect!()

    if plot_needed
        plot!(f, label = "$mod_name $args (constant velocity)")
    end

    println("$mod_name $args (constant velocity): $(Δx*norm(f - f₁))")
    @test Δx*norm(f - f₁) ≈ 0 atol=exp_norm_dev

    advect! = eval(:($mod_name)).generate_solver(f₀, f, args...)
    advect!(v*Δt/Δx)

    if plot_needed
        plot!(f, label = "$mod_name $args")
    end

    println("$mod_name $args: $(Δx*norm(f - f₁))")
    @test Δx*norm(f - f₁) ≈ 0 atol=exp_norm_dev

    if plot_needed
        display(plt)
    end
end

@testset "Test 1D advection solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    v = 1.0

    f₀ = [1.0 + 0.01*sin(2π*i*Δx) for i = 0:100]
    f₁ = [1.0 + 0.01*sin(2π*(i*Δx - v*Δt)) for i = 0:100]

    for mod_name in (:Upwind, :LaxWendroff)
        test_1d_advection_step(mod_name, Δx, Δt,  v, f₀, f₁, 1e-5)
        test_1d_advection_step(mod_name, Δx, Δt, -v, f₁, f₀, 1e-5)
    end

    for riemann_solver in (:Riemann_constant, :Riemann_linear)
        test_1d_advection_step(:Godunov, Δx, Δt,  v, f₀, f₁, 1e-5, riemann_solver)
        test_1d_advection_step(:Godunov, Δx, Δt, -v, f₁, f₀, 1e-5, riemann_solver)
    end

    f₀ = [exp(-((i*Δx - 0.5)/0.15)^2) for i = 0:100]
    f₁ = [exp(-((i*Δx - v*Δt - 0.5)/0.15)^2) for i = 0:100]

    test_1d_advection_step(:Upwind, Δx, Δt,  v, f₀, f₁, 1e-4)
    test_1d_advection_step(:Upwind, Δx, Δt, -v, f₁, f₀, 1e-4)

    test_1d_advection_step(:LaxWendroff, Δx, Δt,  v, f₀, f₁, 1e-5)
    test_1d_advection_step(:LaxWendroff, Δx, Δt, -v, f₁, f₀, 1e-5)

    test_1d_advection_step(:Godunov, Δx, Δt,  v, f₀, f₁, 1e-4, :Riemann_constant)
    test_1d_advection_step(:Godunov, Δx, Δt, -v, f₁, f₀, 1e-4, :Riemann_constant)

    test_1d_advection_step(:Godunov, Δx, Δt,  v, f₀, f₁, 1.1e-4, :Riemann_linear)
    test_1d_advection_step(:Godunov, Δx, Δt, -v, f₁, f₀, 1.1e-4, :Riemann_linear)

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

    test_1d_advection_step(:Upwind, Δx, Δt,  v, f₀, f₁, 1e-18)
    test_1d_advection_step(:Upwind, Δx, Δt, -v, f₀, f₂, 1e-18)

    test_1d_advection_step(:LaxWendroff, Δx, Δt,  v, f₀, f₁, 1.6e-3)
    test_1d_advection_step(:LaxWendroff, Δx, Δt, -v, f₀, f₂, 1.6e-3)

    test_1d_advection_step(:Godunov, Δx, Δt,  v, f₀, f₁, 1e-18, :Riemann_constant)
    test_1d_advection_step(:Godunov, Δx, Δt, -v, f₀, f₂, 1e-18, :Riemann_constant)

    test_1d_advection_step(:Godunov, Δx, Δt,  v, f₀, f₁, 1e-2, :Riemann_linear)
    test_1d_advection_step(:Godunov, Δx, Δt, -v, f₀, f₂, 1e-2, :Riemann_linear)
end
