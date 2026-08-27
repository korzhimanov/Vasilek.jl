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

"""
    test_fdtd_1d_polarization_decoupling(Δx, Δt, cfl)

In 1D the two polarizations are independent: (ey, hz) closes on itself and
so does (ez, hy). Seeding one must leave the other exactly zero forever.
"""
function test_fdtd_1d_polarization_decoupling(Δx, Δt, cfl)
    y = FDTD1D.YeeMesh1D{Float64}(110)
    y.ey[2:102] = [sin(2π*i*Δx) for i = 0:100]
    y.hz[2:101] = [sin(2π*((i+0.5)*Δx-0.5*Δt)) for i = 0:99]
    f = run_fdtd(y, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 10,
                 (t, f) -> zero_current(length(f.ey)))
    println("FDTD1D decoupling y->z: max|ez| = ", maximum(abs, f.ez),
            ", max|hy| = ", maximum(abs, f.hy))
    @test all(iszero, f.ez)
    @test all(iszero, f.hy)

    z = FDTD1D.YeeMesh1D{Float64}(110)
    z.ez[2:102] = [sin(2π*i*Δx) for i = 0:100]
    z.hy[2:101] = [-sin(2π*((i+0.5)*Δx-0.5*Δt)) for i = 0:99]
    f = run_fdtd(z, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 10,
                 (t, f) -> zero_current(length(f.ey)))
    println("FDTD1D decoupling z->y: max|ey| = ", maximum(abs, f.ey),
            ", max|hz| = ", maximum(abs, f.hz))
    @test all(iszero, f.ey)
    @test all(iszero, f.hz)
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

@testset "PML construction" begin
    Δx = 0.01
    Δt = 0.8*Δx

    # As σ_max → 0 the absorbing layer must degenerate into the plain
    # interior update, whose coefficient is the Courant number Δt/Δx.
    # This pins the (Δx, Δt) argument order permanently: swapping them would
    # give Δx/Δt = 1.25 instead of 0.8.
    #
    # σ_max is not free. r₂ = (1-exp(-Δt·σ))/(Δx·σ) loses precision to
    # cancellation when Δt·σ approaches eps, and departs from the linear limit
    # when Δt·σ grows. Measured max relative error against Δt/Δx: 8e-4 at
    # σ_max=1e-8 (cancellation), 5e-7 at σ_max=1.25e-4, 2e-5 at σ_max=1e-6.
    vanishing = FDTD1D.PML(10, 1.25e-4, Δx, Δt)
    @test all(r -> isapprox(r, Δt/Δx; rtol = 1e-5), vanishing.r₂)

    # keyword and positional forms must agree
    positional = FDTD1D.PML(10, 1e3, Δx, Δt)
    keyword = FDTD1D.PML(; N = 10, σ_max = 1e3, Δx = Δx, Δt = Δt)
    @test keyword.r₁ == positional.r₁
    @test keyword.r₂ == positional.r₂

    # and swapping the two must be observable -- otherwise the default
    # argument could stay wrong without any test noticing
    @test FDTD1D.PML(10, 1e3, Δt, Δx).r₂ != positional.r₂
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

    test_fdtd_1d_polarization_decoupling(Δx, Δt, Δt/Δx)

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
