include(joinpath(dirname(@__FILE__), "..","..","src","VlasovSolver","LaxWendroff.jl"))
import .LaxWendroff

include(joinpath(dirname(@__FILE__), "..","..","src","VlasovSolver","Upwind.jl"))
import .Upwind

include(joinpath(dirname(@__FILE__), "..","..","src","VlasovSolver","Godunov.jl"))
import .Godunov

@testset "Test 1D advection solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    v = 1.0
    f₀ = [1.0 + 0.01*sin(2π*i*Δx) for i = 0:1000]
    plt = plot(f₀, label = "initial")
    f₁ = [1.0 + 0.01*sin(2π*(i*Δx - v*Δt)) for i = 0:1000]
    plot!(f₁, label = "expected")

    f = similar(f₀)

    advect! = Upwind.generate_solver(f₀, f, v*Δt/Δx)
    advect!()
    plot!(f, label = "upwind⁺")

    println(Δx*norm(f - f₁))
    @test Δx*norm(f - f₁) ≈ 0 atol=1e-5

    advect! = Upwind.generate_solver(f₁, f, -v*Δt/Δx)
    advect!()
    plot!(f, label = "upwind⁻")

    println(Δx*norm(f - f₀))
    @test Δx*norm(f - f₀) ≈ 0 atol=1e-5

    advect! = Upwind.generate_solver(f₀, f)
    advect!(v*Δt/Δx)
    plot!(f, label = "upwind⁺ c")

    println(Δx*norm(f - f₁))
    @test Δx*norm(f - f₁) ≈ 0 atol=1e-5

    advect! = Upwind.generate_solver(f₁, f)
    advect!(-v*Δt/Δx)
    plot!(f, label = "upwind⁻ c")

    println(Δx*norm(f - f₀))
    @test Δx*norm(f - f₀) ≈ 0 atol=1e-5

    advect! = LaxWendroff.generate_solver(f₀, f)
    advect!(v*Δt/Δx)
    plot!(f, label = "Lax-Wendroff")

    println(Δx*norm(f - f₁))
    @test Δx*norm(f - f₁) ≈ 0 atol=1e-5

    advect! = Godunov.generate_solver(f₀, f, :Riemann_constant)
    advect!(v*Δt/Δx)
    plot!(f, label = "Godunov constant +")

    println(Δx*norm(f - f₁))
    @test Δx*norm(f - f₁) ≈ 0 atol=1e-5

    advect! = Godunov.generate_solver(f₁, f, :Riemann_constant)
    advect!(-v*Δt/Δx)
    plot!(f, label = "Godunov constant -")

    println(Δx*norm(f - f₀))
    @test Δx*norm(f - f₀) ≈ 0 atol=1e-5

    advect! = Godunov.generate_solver(f₀, f, :Riemann_linear)
    advect!(v*Δt/Δx)
    plot!(f, label = "Godunov linear +")

    println(Δx*norm(f - f₁))
    @test Δx*norm(f - f₁) ≈ 0 atol=1e-5
    
    advect! = Godunov.generate_solver(f₁, f, :Riemann_linear)
    advect!(-v*Δt/Δx)
    plot!(f, label = "Godunov linear -")

    println(Δx*norm(f - f₀))
    @test Δx*norm(f - f₀) ≈ 0 atol=1e-5

    xlims!(90, 110)
    display(plt)
end
