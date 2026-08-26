using Vasilek.FDTD1D

const NO_PULSE = (y = (t,x) -> 0.0, z = (t,x) -> 0.0)
no_pml(Δx, Δt) = FDTD1D.PML(0, 1.0, Δx, Δt)

zero_current(n) = (y = zeros(n), z = zeros(n))

function run_fdtd(f₀, cfl, pulse_shape, Δt, Δx, pml, nsteps, j)
    f = deepcopy(f₀)
    advance_fields! = FDTD1D.make_advance_fields(f, cfl, pulse_shape, Δt, Δx, 0, pml)
    for t in 1:nsteps
        advance_fields!(t*Δt, j(t*Δt, f))
    end
    return f
end

function test_fdtd_1d_propagation(Δx, Δt, cfl, f₀, f₁, exp_norm_dev)
    f = run_fdtd(f₀, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 10,
                 (t, f) -> zero_current(length(f.ey)))
    s = sum(@. (f.ey[10:end-1] - f₁.ey[10:end-1])^2)
    println("FDTD1D propagation: $s")
    @test s ≈ 0 atol=exp_norm_dev
end

function test_fdtd_1d_pml(Δx, Δt, cfl, f₀, f₁, exp_norm_dev)
    f = deepcopy(f₀)
    advance_fields! = FDTD1D.make_advance_fields(f, cfl, NO_PULSE, Δt, Δx, 0,
                                                 FDTD1D.PML(10, 1e3, Δx, Δt))
    j = zero_current(length(f.ey))
    for _ in 1:200
        advance_fields!(0, j)
    end
    s = sum(@. (f.ey - f₁.ey)^2)
    println("FDTD1D PML: $s")
    @test s ≈ 0 atol=exp_norm_dev
end

function test_fdtd_1d_generation(Δx, Δt, cfl, f₀, f₁, pulse_shape, exp_norm_dev)
    f = run_fdtd(f₀, cfl, pulse_shape, Δt, Δx, no_pml(Δx, Δt), 100,
                 (t, f) -> zero_current(length(f.ey)))
    s = sum(@. (f.ey - f₁.ey)^2)
    println("FDTD1D generation: $s")
    @test s ≈ 0 atol=exp_norm_dev
end

function test_fdtd_1d_current(Δx, Δt, cfl, f₀, f₁, j, exp_norm_dev)
    f = run_fdtd(f₀, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 100, (t, f) -> j(t))
    s = sum(@. (f.ey - f₁.ey)^2)
    println("FDTD1D current: $s")
    return s
end

@testset "Test 1D FDTD solvers" begin
    Δx = 0.01
    Δt = 0.8*Δx

    f₀ = FDTD1D.YeeMesh1D{Float64}(110)
    f₀.ey[2:102] = [sin(2π*i*Δx) for i = 0:100]
    f₀.hz[2:101] = [sin(2π*((i+0.5)*Δx-0.5*Δt)) for i = 0:99]

    f₁ = FDTD1D.YeeMesh1D{Float64}(110)
    f₁.ey[10:110] = [sin(2π*(i*Δx)) for i = 0:100]
    f₁.hz[10:109] = [sin(2π*((i+0.5)*Δx-0.5*Δt)) for i = 0:99]

    test_fdtd_1d_propagation(Δx, Δt, Δt/Δx, f₀, f₁, 1e-3)

    f₁ = FDTD1D.YeeMesh1D{Float64}(110)
    test_fdtd_1d_pml(Δx, Δt, Δt/Δx, f₀, f₁, 1e-5)

    f₀ = FDTD1D.YeeMesh1D{Float64}(100)
    f₁ = FDTD1D.YeeMesh1D{Float64}(100)
    f₁.ey[3:81] = [-sin(2π*(i*Δx-100*Δt)) for i = 1:79]
    f₁.hz[3:81] = [-sin(2π*((i+0.5)*Δx-100.5*Δt)) for i = 1:79]

    pulse_shape = (y = (t,x) -> sin(2π*(x-t)), z = (t,x) -> 0.0)
    test_fdtd_1d_generation(Δx, Δt, Δt/Δx, f₀, f₁, pulse_shape, 0.02)

    f₀ = FDTD1D.YeeMesh1D{Float64}(200)
    f₁ = FDTD1D.YeeMesh1D{Float64}(200)
    f₁.ey[100:178] = [-sin(2π*(i*Δx-100*Δt)) for i = 0:78]
    f₁.ey[22:100] = [-sin(2π*(i*Δx-100*Δt)) for i = 78:-1:0]
    f₁.hz[3:81] = [-sin(2π*((i+0.5)*Δx-100.5*Δt)) for i = 1:79]

    function j(t)
        jy = zeros(length(f₀.ey))
        jy[end÷2] = π/2*sin(2π*t)
        return (y = jy, z = zeros(length(f₀.ey)))
    end

    # The assertion here was commented out when this test was written, so the
    # function only ever printed. It passes -- a working check had been left
    # disabled.
    s = test_fdtd_1d_current(Δx, Δt, Δt/Δx, f₀, f₁, j, 0.1)
    @test s ≈ 0 atol=0.1
end
