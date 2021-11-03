include(joinpath("..","..","src","MaxwellSolver","FDTD1D.jl"))
import .FDTD1D

function test_fdtd_1d_propagation(Δx, Δt, cfl, f₀, f₁, exp_norm_dev, args...; plot_needed=false, kwargs...)
    if plot_needed
        plt = plot(f₀.ey, label = "initial")
        plot!(f₁.ey, label = "expected")
    end

    shape(t,x) = 0.0

    pulse_shape = Dict(
        "y"=>shape,
        "z"=>shape
    )

    f = deepcopy(f₀)

    advance_fields! = FDTD1D.make_advance_fields(f, cfl, pulse_shape, Δt, Δx, 0)
    
    for t in 1:10
        advance_fields!(0, :Nothing)
    end

    if plot_needed
        plot!(collect(11:f.N-9), f.ey[11:end-10], label = "calculated")
    end

    s = sum(@. (f.ey[11:end-10] - f₁.ey[11:end-10])^2)
    println("FDTD1D propagation $args $(values(values(kwargs))): $(s)")
    @test s ≈ 0 atol=exp_norm_dev

    if plot_needed
        display(plt)
    end
end

function test_fdtd_1d_generation(Δx, Δt, cfl, f₀, f₁, pulse_shape, exp_norm_dev, args...; plot_needed=false, kwargs...)
    if plot_needed
        plt = plot(f₀.ey, label = "initial")
        plot!(f₁.ey, label = "expected")
    end

    f = deepcopy(f₀)

    advance_fields! = FDTD1D.make_advance_fields(f, cfl, pulse_shape, Δt, Δx, 0)
    
    for t in 1:100
        advance_fields!(t*Δt, :Nothing)
    end

    if plot_needed
        plot!(f.ey, label = "calculated")
    end

    s = sum(@. (f.ey - f₁.ey)^2)
    println("FDTD1D generation $args $(values(values(kwargs))): $(s)")
    @test s ≈ 0 atol=exp_norm_dev

    if plot_needed
        display(plt)
    end
end

@testset "Test 1D FDTD solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx
    
    f₀ = FDTD1D.YeeMesh1D{Float64}(100)
    f₀.ey[:] = [sin(2π*i*Δx) for i = 0:100]
    f₀.hz[:] = [sin(2π*((i+0.5)*Δx-0.5*Δt)) for i = 0:99]

    f₁ = FDTD1D.YeeMesh1D{Float64}(100)
    f₁.ey[:] = [sin(2π*(i*Δx-10*Δt)) for i = 0:100]
    f₁.hz[:] = [sin(2π*((i+0.5)*Δx-10.5*Δt)) for i = 0:99]

    test_fdtd_1d_propagation(Δx, Δt, Δt/Δx, f₀, f₁, 1e-7)

    f₀ = FDTD1D.YeeMesh1D{Float64}(100)
    
    f₁ = FDTD1D.YeeMesh1D{Float64}(100)
    f₁.ey[3:81] = [-sin(2π*(i*Δx-100*Δt)) for i = 1:79]
    f₁.hz[3:81] = [-sin(2π*((i+0.5)*Δx-100.5*Δt)) for i = 1:79]
    
    shape_y(t,x) = sin(2π*(x-t))
    shape_z(t,x) = 0.0

    pulse_shape = Dict(
        "y"=>shape_y,
        "z"=>shape_z
    )

    test_fdtd_1d_generation(Δx, Δt, Δt/Δx, f₀, f₁, pulse_shape, 0.02)
    # plot!(f.ey, label = "calculated")
end
