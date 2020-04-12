include(joinpath("..","..","src","VlasovSolver","LaxWendroff.jl"))
import .LaxWendroff

@testset "Test 1D advection solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    v = 1
    f₀ = [1.0 + 0.01*sin(2π*i*Δx) for i = 0:1000]
    f = similar(f₀)
    advect! = LaxWendroff.generate_solver(f₀, f, v*Δt/Δx)
    advect!()

    plt = plot(f₀)
    plot!(f)
    plot!([1.0 + 0.01*sin(2π*(i*Δx - v*Δt)) for i = 0:1000])
    xlims!(90, 110)
    display(plt)

    @test Δx*norm(f - [1.0 + 0.01*sin(2π*i*Δx - v*Δt) for i = 0:1000]) ≈ 0 atol=1e-4
end
