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

    for mod_name in (:LaxWendroff, :Upwind)
        test_1d_advection_step(mod_name, Δx, Δt,  v, f₀, f₁, 1e-5)
        test_1d_advection_step(mod_name, Δx, Δt, -v, f₁, f₀, 1e-5)
    end

    for riemann_solver in (:Riemann_constant, :Riemann_linear)
        test_1d_advection_step(:Godunov, Δx, Δt,  v, f₀, f₁, 1e-5, riemann_solver)
        test_1d_advection_step(:Godunov, Δx, Δt, -v, f₁, f₀, 1e-5, riemann_solver)
    end
end
