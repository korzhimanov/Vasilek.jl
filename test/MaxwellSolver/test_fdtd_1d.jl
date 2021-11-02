include(joinpath("..","..","src","MaxwellSolver","FDTD1D.jl"))
import .FDTD1D

@testset "Test 1D FDTD solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    
    f₀ = FDTD1D.YeeMesh1D{Float64}(100)
    f₀.ey[:] = [sin(2π*i*Δx) for i = 0:100]
    f₀.hz[:] = [sin(2π*((i+0.5)*Δx-0.5*Δt)) for i = 0:99]

    f₁ = FDTD1D.YeeMesh1D{Float64}(100)
    f₁.ey[:] = [sin(2π*(i*Δx-10*Δt)) for i = 0:100]
    f₁.hz[:] = [sin(2π*((i+0.5)*Δx-10.5*Δt)) for i = 0:99]

    shape(t,x) = 0.0

    pulse_shape = Dict(
        "y"=>shape,
        "z"=>shape
    )

    plt = plot(f₀.ey, label = "initial")
    plot!(f₁.ey, label = "expected")

    f = deepcopy(f₀)

    advance_fields! = FDTD1D.make_advance_fields(f, Δt/Δx, pulse_shape, Δt, Δx, 0)

    for t in 1:10
        advance_fields!(0, :Nothing)
    end

    plot!([11:90], f.ey[11:90], label = "calculated")
    plot!(f.ey, label = "calculated")
    @test sum(@. (f.ey[11:90] - f₁.ey[11:90])^2) ≈ 0 atol=1e-7
end
