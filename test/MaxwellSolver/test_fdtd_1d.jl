using Vasilek.FDTD1D

const NO_PULSE = (y = (t,x) -> 0.0, z = (t,x) -> 0.0)
no_pml(Δx, Δt) = FDTD1D.PML(0, 1.0, Δx, Δt)

zero_current(n) = (y = zeros(n), z = zeros(n))

function run_fdtd(f₀, cfl, pulse_shape, Δt, Δx, pml, nsteps, j)
    f = deepcopy(f₀)
    op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = pulse_shape, pml)
    for t in 1:nsteps
        FDTD1D.advance!(f, op, t*Δt, j(t*Δt, f))
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
    op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = FDTD1D.PML(10, 1e3, Δx, Δt))
    j = zero_current(length(f.ey))
    for _ in 1:200
        FDTD1D.advance!(f, op, 0, j)     # an integer time, as a source may be handed
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
    # r₂ = (1-exp(-Δt·σ))/(Δx·σ) is taken through `expm1`, so it no longer
    # cancels as Δt·σ → 0: written as `1 - exp`, σ_max = 1e-9 gave r₂ between
    # 0.4974 and 0.50002 of a 0.5, and σ_max = 0 gave 0/0.
    vanishing = FDTD1D.PML(10, 1.25e-4, Δx, Δt)
    @test all(r -> isapprox(r, Δt/Δx; rtol = 1e-5), vanishing.r₂)
    for σ_max in (1e-9, 1e-12, 0.0)
        tiny = FDTD1D.PML(10, σ_max, Δx, Δt)
        # the exact value departs from Δt/Δx by Δt·σ/2 ≤ 4e-12 at 1e-9
        @test all(r -> isapprox(r, Δt/Δx; rtol = 1e-10), tiny.r₂)
    end
    @test_throws ArgumentError FDTD1D.PML(10, -1.0, Δx, Δt)

    # keyword and positional forms must agree
    positional = FDTD1D.PML(10, 1e3, Δx, Δt)
    keyword = FDTD1D.PML(; N = 10, σ_max = 1e3, Δx = Δx, Δt = Δt)
    @test keyword.r₁ == positional.r₁
    @test keyword.r₂ == positional.r₂

    # and swapping the two must be observable -- otherwise the default
    # argument could stay wrong without any test noticing
    @test FDTD1D.PML(10, 1e3, Δt, Δx).r₂ != positional.r₂

    # An integer σ_max is a number like any other. It was a MethodError while
    # the type took its float parameter from σ_max alone.
    whole = FDTD1D.PML(10, 1000, Δx, Δt)
    @test whole isa FDTD1D.PML{Float64}
    @test whole.σ_max === 1e3
    @test whole.r₁ == positional.r₁
    @test whole.r₂ == positional.r₂
    @test FDTD1D.PML(; N = 10, σ_max = 1000, Δx, Δt).r₂ == positional.r₂
    @test FDTD1D.PML(0, 1, Δx, Δt) isa FDTD1D.PML{Float64}

    # A wider type than Float64 computes its depth profile in that type, so
    # the layer is as accurate as its type: (i/2N)^3 was always Float64.
    wider = FDTD1D.PML(3, big"1e3", big"0.01", big"0.008")
    σ = [big"1e3"*(BigFloat(i)/6)^3 for i = 1:6]
    @test wider isa FDTD1D.PML{BigFloat}
    @test maximum(abs, wider.r₁ .- exp.(-big"0.008" .* σ)) < 1e-70
end

@testset "Yee1D checks its arguments" begin
    Δx = 0.01; Δt = 0.8*Δx
    m = FDTD1D.YeeMesh1D{Float64}(50)
    pml = FDTD1D.PML(; N = 10, σ_max = 1e3, Δx = Δx, Δt = Δt)
    # the interior would run at 0.9, the layer at 0.8
    @test_throws ArgumentError FDTD1D.Yee1D(; Δx, Δt, cfl = 0.9, source = NO_PULSE, pml)
    # and so does the positional constructor, which is the one that checks
    @test_throws ArgumentError FDTD1D.Yee1D(0.9, Δx, Δt, 0.0, pml, NO_PULSE)
    @test FDTD1D.Yee1D(0.8, Δx, Δt, 0.0, pml, NO_PULSE).cfl === 0.8
    # The step refuses a mesh or a current it cannot take before it writes
    # anything. Every field is seeded and the source is on, so a check that ran
    # after the injection or after either polarisation would show.
    pulse = (y = (t, x) -> 1.0, z = (t, x) -> 1.0)
    function refuses_untouched(m, op, j)
        m.ey[26] = m.ez[26] = m.hy[25] = m.hz[25] = 1.0
        before = deepcopy(m)
        thrown = try
            FDTD1D.advance!(m, op, 0.0, j); nothing
        catch e
            e
        end
        @test m.ey == before.ey && m.ez == before.ez &&
              m.hy == before.hy && m.hz == before.hz
        return thrown
    end
    # two layers of 30 cells do not fit in 50: the operator is built without a
    # mesh, so it is the step that refuses
    wide = FDTD1D.Yee1D(; Δx, Δt, source = pulse,
                        pml = FDTD1D.PML(; N = 30, σ_max = 1e3, Δx = Δx, Δt = Δt))
    @test refuses_untouched(m, wide, zero_current(51)) isa ArgumentError
    @test FDTD1D.advance!(FDTD1D.YeeMesh1D{Float64}(62), wide, 0.0, zero_current(63)) isa
          FDTD1D.YeeMesh1D{Float64}                  # 2·30 + 2 cells is enough
    # a current on fewer (or more) nodes than the mesh has
    driven = FDTD1D.Yee1D(; Δx, Δt, source = pulse, pml)
    @test refuses_untouched(m, driven, zero_current(26)) isa DimensionMismatch
    @test refuses_untouched(m, driven, (y = zeros(51), z = zeros(26))) isa DimensionMismatch
    @test refuses_untouched(m, driven, zero_current(52)) isa DimensionMismatch
    # the axes, not the length: a vector of 51 entries indexed from 0 passed
    # the length check and threw on `j[51]` with the source already injected
    shifted = (y = Base.IdentityUnitRange(0:50), z = Base.IdentityUnitRange(0:50))
    @test refuses_untouched(m, driven, shifted) isa DimensionMismatch
    fill!(m.ey, 0); fill!(m.ez, 0); fill!(m.hy, 0); fill!(m.hz, 0)

    # the defaults: the Courant number of the step, no offset, the layer the
    # solver has always defaulted to
    op = FDTD1D.Yee1D(; Δx, Δt, source = NO_PULSE)
    @test op isa FDTD1D.Yee1D{Float64}
    @test op.cfl == Δt/Δx
    @test op.x_min === 0.0
    @test op.pml.N == 10 && op.pml.σ_max == 1e3
    @test op.pml.r₂ == pml.r₂
    # an integer x_min, as every caller passed before, is a coordinate
    @test FDTD1D.Yee1D(; Δx, Δt, source = NO_PULSE, x_min = 0).x_min === 0.0
    # cfl is kept as given, not recomputed: Δt/Δx need not round back to it
    @test FDTD1D.Yee1D(; Δx, Δt, cfl = 0.8, source = NO_PULSE).cfl === 0.8

    # The step works in the type of Δx and Δt. A Float64 literal cfl, the
    # default layer's Float64 σ_max and a Float64 layer passed in are rounded
    # to it; cfl is checked there, where Δt/Δx is one ulp off 0.8.
    Δx32 = 0.01f0; Δt32 = 0.008f0
    @test Δt32/Δx32 != 0.8f0
    for kw in ((;), (cfl = 0.8,), (cfl = 0.8f0,), (pml = pml,))
        op32 = FDTD1D.Yee1D(; Δx = Δx32, Δt = Δt32, source = NO_PULSE, kw...)
        @test op32 isa FDTD1D.Yee1D{Float32}
        @test op32.pml isa FDTD1D.PML{Float32}
    end
    @test FDTD1D.Yee1D(; Δx = Δx32, Δt = Δt32, cfl = 0.8, source = NO_PULSE).cfl === 0.8f0
    @test_throws ArgumentError FDTD1D.Yee1D(; Δx = Δx32, Δt = Δt32, cfl = 0.8001,
                                             source = NO_PULSE)
    # the layer is rounded once, from the coefficients computed in Float64
    @test FDTD1D.Yee1D(; Δx = Δx32, Δt = Δt32, source = NO_PULSE).pml.r₁ ==
          Float32.(FDTD1D.PML(; N = 10, σ_max = 1e3, Δx = Δx32, Δt = Δt32).r₁)
    # x_min is a coordinate and keeps its type: the source sees x_min + Δx in
    # Float64, as it did before the step had a type
    xs = Float64[]
    probe = (y = (t, x) -> (push!(xs, x); 0.0), z = (t, x) -> 0.0)
    op32 = FDTD1D.Yee1D(; Δx = Δx32, Δt = Δt32, x_min = -5*2π, source = probe)
    @test op32.x_min === -5*2π
    FDTD1D.advance!(FDTD1D.YeeMesh1D{Float32}(50), op32, 0.0f0,
                    (y = zeros(Float32, 51), z = zeros(Float32, 51)))
    @test xs == [-5*2π + Δx32, -5*2π + 1.5*Δx32]

    # Both injections call the source with a time in the step's type: the
    # magnetic one was `t + 0.5*Δt`, which the literal widened to Float64.
    ts = DataType[]
    timed = (y = (t, x) -> (push!(ts, typeof(t)); 0.0f0), z = (t, x) -> 0.0f0)
    FDTD1D.advance!(FDTD1D.YeeMesh1D{Float32}(50),
                    FDTD1D.Yee1D(; Δx = Δx32, Δt = Δt32, source = timed), 0.0f0,
                    (y = zeros(Float32, 51), z = zeros(Float32, 51)))
    @test ts == [Float32, Float32]

    # Steps that are not finite and positive gave Inf, NaN or a scheme running
    # backwards in time; a layer built for Δt and Δx swapped, an edge that
    # reflects. Both used to be accepted.
    for (h, τ) in ((0.0, 0.008), (-0.01, -0.008), (Inf, Inf), (0.01, NaN))
        @test_throws ArgumentError FDTD1D.Yee1D(; Δx = h, Δt = τ, source = NO_PULSE,
                                                 pml = no_pml(0.01, 0.008))
    end
    @test_throws ArgumentError FDTD1D.Yee1D(; Δx, Δt, source = NO_PULSE,
                                             pml = FDTD1D.PML(10, 1e3, Δt, Δx))
    # with σ_max = 0 the layer's r₂ is Δt/Δx itself, and swapped it is not
    @test_throws ArgumentError FDTD1D.Yee1D(; Δx, Δt, source = NO_PULSE,
                                             pml = FDTD1D.PML(10, 0.0, Δt, Δx))
    # The operator copies the layer: it is an immutable value, and shared
    # coefficients would change under it with the caller's.
    @test op.pml.r₁ !== pml.r₁ && op.pml.r₁ == pml.r₁
    shared = FDTD1D.Yee1D(; Δx, Δt, source = NO_PULSE, pml)
    @test shared.pml.r₁ !== pml.r₁ && shared.pml.r₂ !== pml.r₂

    # no scratch, whatever generic code passes along
    @test workspace(op) === nothing
    @test workspace(op, 50) === nothing

    # a step returns the mesh it advanced, as `advect!` and `solve!` return
    # their destinations, and infers
    j = zero_current(51)
    @test FDTD1D.advance!(m, op, 0.0, j) === m
    @test (@inferred FDTD1D.advance!(m, op, 0.0, j)) === m
end

@testset "Both polarisations are launched rightwards" begin
    # (ez, hy) obeys the (ey, hz) equations with hy = -hz, so the same pulse in
    # `z` must produce the same ez as in `y` it produces ey, bit for bit. The z
    # source used to inject hy with the sign of hz, which sent all but 0.8% of
    # the pulse left into the absorbing layer.
    Δx = 0.05; Δt = 0.5*Δx; N = 400
    pulse(t, x) = t < 2 ? sin(2π*(x - t))*sin(π*t/2)^2 : 0.0
    nothing_(t, x) = 0.0
    pml = FDTD1D.PML(; N = 20, σ_max = 1e3, Δx = Δx, Δt = Δt)
    j(t, f) = zero_current(length(f.ey))
    fy = run_fdtd(FDTD1D.YeeMesh1D{Float64}(N), Δt/Δx, (y = pulse, z = nothing_),
                  Δt, Δx, pml, 200, j)
    fz = run_fdtd(FDTD1D.YeeMesh1D{Float64}(N), Δt/Δx, (y = nothing_, z = pulse),
                  Δt, Δx, pml, 200, j)
    println("FDTD1D one-way source: Σey² = ", sum(abs2, fy.ey), ", Σez² = ", sum(abs2, fz.ez))
    @test sum(abs2, fy.ey) > 1          # the pulse is on the grid, not in the layer
    @test fz.ez == fy.ey
    @test fz.hy == -fy.hz
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

"A right-going Gaussian, seeded consistently for the given Courant number."
function gaussian_mesh(N, Δx, cfl; centre = 60, width = 3.0, halfwidth = 20)
    m = FDTD1D.YeeMesh1D{Float64}(N)
    for i = 0:2*halfwidth
        m.ey[centre - halfwidth + i] = exp(-((i - halfwidth)/width)^2)
        m.hz[centre - halfwidth + i] = exp(-((i + 0.5 - halfwidth - 0.5*cfl)/width)^2)
    end
    return m
end

@testset "FDTD at the magic time step is exact" begin
    # At cfl = 1 the 1D Yee scheme has no numerical dispersion at all: the
    # update reduces to a one-cell shift per step, and a pulse translates
    # exactly. This is the strongest statement available about the interior
    # update, and it is a near-bit-level one -- measured 1.1e-19 after 20 steps
    # against a plain circshift, where a single mis-signed or mis-indexed term
    # would leave an O(1) residue.
    #
    # cfl = 0.8 is run alongside to show the test has teeth: the same
    # comparison there is off by 0.83, which is the physical dispersion the
    # existing propagation test tolerates.
    Δx = 0.01; N = 200
    for (cfl, exact) in ((1.0, true), (0.8, false))
        Δt = cfl*Δx
        m = gaussian_mesh(N, Δx, cfl)
        before = copy(m.ey)
        f = run_fdtd(m, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 20,
                     (t, f) -> zero_current(length(f.ey)))
        dev = maximum(abs, f.ey[80:120] .- circshift(before, 20)[80:120])
        println("  cfl = ", cfl, ": max|ey − shift(ey₀, 20)| = ", dev)
        if exact
            @test dev < 1e-15
        else
            @test dev > 0.1        # dispersion is real, and the test can see it
        end
    end
end

@testset "FDTD Courant stability limit" begin
    # cfl ≤ 1 is the 1D stability condition, and nothing asserted it. Seeded
    # from a single nonzero cell so every mode including the grid-scale one is
    # excited. Measured after 2000 steps: 0.409 at cfl = 0.99, exactly 1.0 at
    # cfl = 1, and NaN at 1.02 and 1.2.
    Δx = 0.01; N = 200
    for cfl in (0.99, 1.0)
        Δt = cfl*Δx
        m = FDTD1D.YeeMesh1D{Float64}(N); m.ey[100] = 1.0
        f = run_fdtd(m, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 2000,
                     (t, f) -> zero_current(length(f.ey)))
        peak = maximum(abs, f.ey)
        println("  cfl = ", rpad(cfl, 5), " after 2000 steps: max|ey| = ", peak)
        @test all(isfinite, f.ey)
        @test peak ≤ 1.0 + 1e-12
    end
    for cfl in (1.02, 1.2)
        Δt = cfl*Δx
        m = FDTD1D.YeeMesh1D{Float64}(N); m.ey[100] = 1.0
        f = run_fdtd(m, cfl, NO_PULSE, Δt, Δx, no_pml(Δx, Δt), 2000,
                     (t, f) -> zero_current(length(f.ey)))
        println("  cfl = ", rpad(cfl, 5), " after 2000 steps: max|ey| = ", maximum(abs, f.ey))
        @test !all(isfinite, f.ey)
    end
end

@testset "YeeMesh1D shape" begin
    # The staggering: N+1 electric nodes, N magnetic cells between them. Nothing
    # asserted it, and every index expression in the module depends on it.
    for T in (Float64, Float32)
        m = FDTD1D.YeeMesh1D{T}(7)
        @test typeof(m) === FDTD1D.YeeMesh1D{T}    # one parameter: N is an Int
        @test length(m.ey) == length(m.ez) == 8
        @test length(m.hy) == length(m.hz) == 7
        @test m.N == 7
        @test FDTD1D.YeeMesh1D{T}(Int32(7)).N === 7
        @test eltype(m.ey) === eltype(m.hz) === T
        @test all(iszero, m.ey) && all(iszero, m.ez)
        @test all(iszero, m.hy) && all(iszero, m.hz)
    end
end

@testset "The Yee leapfrog conserves its staggered energy" begin
    # `E` sits at integer steps and `H` at half-integer ones, so the conserved
    # quadratic form is staggered in time too:
    #
    #     U = ‖E^{n+1}‖² + ⟨H^{n+1/2}, H^{n+3/2}⟩
    #
    # which the update makes exact. Writing `E^{n+1} = E^n − A H^{n+1/2}` and
    # `H^{n+3/2} = H^{n+1/2} + Aᵀ E^{n+1}` -- the discrete curls are adjoint,
    # the boundary terms vanishing because `ey[1]` and `ey[end]` are frozen --
    # gives `U^{n+1} = U^n` identically, and the equivalent form
    # `⟨E^n, E^{n+1}⟩ + ‖H^{n+1/2}‖²`.
    #
    # The naive `‖E‖² + ‖H‖²` at a single instant is *not* conserved and must
    # not be used: measured relative range over 2000 steps is 0.12 at cfl = 0.5,
    # 0.20 at 0.8 and 0.24 at 0.99, against 1e-15 for the staggered form. Both
    # staggered forms agree to 9e-15.
    Δx = 0.01; N = 400
    for cfl in (0.5, 0.8, 0.99)
        Δt = cfl*Δx
        m = FDTD1D.YeeMesh1D{Float64}(N)
        for i = 0:N
            m.ey[i+1] = exp(-((i - 200)/20)^2)*sin(2π*i/25)
        end
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = no_pml(Δx, Δt))
        j = zero_current(N + 1)

        staggered = Float64[]; equivalent = Float64[]; naive = Float64[]
        for s = 1:2000
            e_pre = copy(m.ey); h_pre = copy(m.hz)
            FDTD1D.advance!(m, op, s*Δt, j)
            push!(naive,      sum(abs2, e_pre) + sum(abs2, h_pre))
            push!(staggered,  sum(abs2, m.ey) + sum(h_pre .* m.hz))
            push!(equivalent, sum(e_pre .* m.ey) + sum(abs2, h_pre))
        end
        spread(u) = (maximum(u) - minimum(u))/abs(u[1])
        println("  cfl = ", rpad(cfl, 5), " staggered energy drift = ",
                rpad(round(spread(staggered); sigdigits = 3), 11),
                " naive = ", round(spread(naive); sigdigits = 3))
        @test spread(staggered) < 1e-12
        @test spread(equivalent) < 1e-12
        @test maximum(abs, staggered .- equivalent) < 1e-12
        @test spread(naive) > 0.05          # and the naive form genuinely is not
    end
end

@testset "PML absorbs, equally at both ends" begin
    # The existing PML test checks the field ends up near zero. That passes for
    # a layer that absorbs badly but symmetrically, and for one that reflects
    # into a mode the tolerance happens not to see. Measured here instead: the
    # reflection coefficient, and the left/right asymmetry.
    #
    # The index arithmetic differs between the two ends --
    # `r₁[1+2*(N-i+1)]` against `r₁[1+2*(i-Nx+N-1)]` -- so an off-by-one in one
    # of them is invisible to a test that only looks at one direction.
    #
    # Measured: R = 6.25e-9 rightgoing, 6.26e-9 leftgoing, asymmetry 1.0010.
    # The same run without a PML leaves 0.999 of the pulse bouncing around.
    Δx = 0.01; cfl = 0.8; Δt = cfl*Δx; N = 300; NP = 10

    function residual(direction, pml)
        m = FDTD1D.YeeMesh1D{Float64}(N)
        for i = 0:N
            m.ey[i+1] = exp(-((i - 150)/12)^2)
            i + 1 ≤ N && (m.hz[i+1] = direction*exp(-((i + 0.5 - 150 - 0.5*cfl)/12)^2))
        end
        incident = maximum(abs, m.ey)
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml)
        j = zero_current(N + 1)
        for s = 1:600
            FDTD1D.advance!(m, op, s*Δt, j)
        end
        return maximum(abs, m.ey[NP+2:N-NP])/incident
    end

    absorbing = FDTD1D.PML(NP, 1e3, Δx, Δt)
    right = residual(+1.0, absorbing)
    left  = residual(-1.0, absorbing)
    println("  reflection: rightgoing ", right, ", leftgoing ", left,
            ", asymmetry ", max(left, right)/min(left, right))
    @test right < 1e-7
    @test left < 1e-7
    @test max(left, right)/min(left, right) < 1.05

    # Without the layer the pulse is still there, so the numbers above are the
    # PML working rather than the pulse having left the grid.
    @test residual(+1.0, no_pml(Δx, Δt)) > 0.5
end

@testset "A vanishing PML leaves the interior alone" begin
    # As σ_max → 0 the layer degenerates into the plain update, so the interior
    # must evolve as if there were no layer at all. The existing PML test checks
    # this on the coefficient `r₂`; this checks it on the field, which is what
    # actually matters. Measured max|Δ| over 50 steps: 2.4e-17 at σ_max = 1.25e-4.
    Δx = 0.01; cfl = 0.8; Δt = cfl*Δx; N = 200
    function run(pml)
        m = FDTD1D.YeeMesh1D{Float64}(N)
        for i = 0:N
            m.ey[i+1] = exp(-((i - 100)/12)^2)
        end
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml)
        j = zero_current(N + 1)
        for s = 1:50
            FDTD1D.advance!(m, op, s*Δt, j)
        end
        return m
    end
    bare = run(no_pml(Δx, Δt))
    faint = run(FDTD1D.PML(10, 1.25e-4, Δx, Δt))
    dev = maximum(abs, faint.ey[15:N-14] .- bare.ey[15:N-14])
    println("  σ_max = 1.25e-4 vs no layer: max|Δ| in the interior = ", dev)
    @test dev < 1e-12
end

@testset "The end nodes are PEC boundaries" begin
    # Tangential E vanishes at a perfect conductor, and this solver imposes that
    # by never writing to either end node: the interior loop runs
    # `pml.N+2 : Nx-pml.N` and the two absorbing-layer loops stop short of both
    # ends, so `ey[1]` and `ey[end]` hold the zero `YeeMesh1D` gives them.
    #
    # They are not dead storage -- the `hz` update (`_update_h!`) reads
    # `ey[end]` -- so this is the boundary condition rather than an accident of
    # the loop bounds, and it is what makes the staggered energy above exactly
    # conserved: the discrete curls are adjoint only because the boundary terms
    # vanish.
    #
    # The current used to be added over `1:Nx`, which drove node 1 while node
    # `Nx+1` was left alone. Injecting a current into a perfect conductor is
    # meaningless, and doing it at one end only broke the symmetry the energy
    # identity depends on. It is now applied over `2:Nx`, exactly the set of
    # dynamic nodes.
    Δx = 0.01; cfl = 0.8; Δt = cfl*Δx; N = 20

    @testset "a current drives the interior only, PML N = $NP" for NP in (0, 5)
        m = FDTD1D.YeeMesh1D{Float64}(N)
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = FDTD1D.PML(NP, 1e3, Δx, Δt))
        FDTD1D.advance!(m, op, 0.0, (y = fill(1.0, N+1), z = fill(1.0, N+1)))
        @test length(m.ey) == N + 1
        @test m.ey[1] == 0.0                 # PEC: not driven
        @test m.ey[N+1] == 0.0               # PEC: not driven
        @test all(m.ey[2:N] .== 1.0)         # every interior node is
        @test m.ez[1] == 0.0 && m.ez[N+1] == 0.0
        @test all(m.ez[2:N] .== 1.0)
    end

    @testset "the ends stay zero through a real run, PML N = $NP" for NP in (0, 10)
        m = FDTD1D.YeeMesh1D{Float64}(200)
        for i = 1:120
            m.ey[i+40] = exp(-((i - 60)/12)^2)     # seeded in the interior
        end
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = FDTD1D.PML(NP, 1e3, Δx, Δt))
        j = (y = fill(1e-3, 201), z = fill(1e-3, 201))
        for s = 1:300
            FDTD1D.advance!(m, op, s*Δt, j)
        end
        @test m.ey[1] == 0.0
        @test m.ey[end] == 0.0
        @test m.ez[1] == 0.0
        @test m.ez[end] == 0.0
        @test any(!iszero, m.ey)             # and the run did something
    end

    @testset "a wave reflects with inverted sign, PML N = 0" begin
        # The observable consequence of PEC, and the reason it matters: a pulse
        # hitting the wall comes back inverted. An open or absorbing end would
        # not do this, so it distinguishes the boundary condition rather than
        # merely observing that two array slots stay zero.
        N = 400
        m = FDTD1D.YeeMesh1D{Float64}(N)
        for i = 0:N
            m.ey[i+1] = exp(-((i - 300)/12)^2)
            i + 1 ≤ N && (m.hz[i+1] = exp(-((i + 0.5 - 300 - 0.5*cfl)/12)^2))
        end
        incident = maximum(m.ey)
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = no_pml(Δx, Δt))
        j = zero_current(N + 1)
        for s = 1:250                        # out to the wall and part way back
            FDTD1D.advance!(m, op, s*Δt, j)
        end
        println("  PEC reflection: incident ", incident, ", reflected ", minimum(m.ey),
                ", ratio ", minimum(m.ey)/incident)
        @test minimum(m.ey)/incident < -0.9  # inverted, and almost lossless
        @test maximum(m.ey) < 0.1*incident   # nothing of the original sign left
    end
end

@testset "The Yee scheme's numerical dispersion relation" begin
    # The 1D Yee update satisfies, exactly,
    #
    #     sin(ωΔt/2) = cfl · sin(kΔx/2)
    #
    # and nothing measured it. The suite asserts the magic-step case (`cfl = 1`,
    # no dispersion at all) and, at `cfl = 0.8`, only that the deviation from a
    # pure translation *exceeds* 0.1 -- the error is bounded from below and not
    # from above, so a scheme that was wrong but dispersive would pass. This
    # closes it from the other side, against a closed form rather than a
    # tolerance.
    #
    # A PEC standing mode `sin(kx)` with `k = mπ/L` is the natural probe here:
    # both end nodes are held at zero by the boundary condition, so the mode is
    # an exact eigenfunction of the discrete operator and oscillates as
    # `cos(ωt)`.
    #
    # **Measured through the mode's own projection `Σ ey·sin(kx)`, not a point
    # sample.** A single probe at `L/4` reads `sin(mπ/4)`, which is exactly zero
    # for every `m` divisible by 4 -- at `m = 20` the "signal" is round-off, and
    # the fitted frequency came out 4.2x too high. The projection has no nodes.
    #
    # Measured relative departure from the closed form: 1.4e-16 to 1.3e-4 over
    # `cfl ∈ {0.5, 0.9, 1.0}` and `kΔx` from 0.016 to 1.41.
    Δx = 0.01
    N = 200
    L = N*Δx
    worst = 0.0

    for cfl in (0.5, 0.9, 1.0), m in (1, 5, 20, 50, 90)
        Δt = cfl*Δx
        k = m*π/L
        mesh = FDTD1D.YeeMesh1D{Float64}(N)
        shape = [sin(k*i*Δx) for i = 0:N]
        mesh.ey .= shape
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = no_pml(Δx, Δt))
        j = zero_current(N + 1)
        amplitude = Float64[]
        for s = 1:4000
            FDTD1D.advance!(mesh, op, s*Δt, j)
            push!(amplitude, sum(mesh.ey .* shape))
        end

        ups = [s for s in 2:length(amplitude)
               if amplitude[s-1] < 0 ≤ amplitude[s]]
        ω = 2π/((ups[end] - ups[1])/(length(ups) - 1)*Δt)
        ω_yee = 2/Δt*asin(min(1.0, cfl*sin(k*Δx/2)))
        dev = abs(ω - ω_yee)/ω_yee
        worst = max(worst, dev)
        @test dev < 1e-3

        # At the magic time step the relation collapses to ω = ck exactly:
        # `asin(sin(kΔx/2))·2/Δt = k·Δx/Δt = k`. Measured departure from the
        # *ideal* dispersion 0.0 at every mode but the last, where the mode is
        # near the grid scale and the crossing count quantises: 1.1e-4.
        if cfl == 1.0
            @test abs(ω_yee - k)/k < 1e-14
        end
    end
    println("  worst departure from sin(ωΔt/2) = cfl·sin(kΔx/2): ", worst)

    # The physical content, which the closed form makes quantitative: short
    # waves travel slow. Measured phase velocity at kΔx = 1.41 (λ ≈ 4.4Δx):
    # 0.936c at cfl = 0.5, 0.981c at cfl = 0.9, and exactly c at cfl = 1.
    for (cfl, expected) in ((0.5, 0.93575), (0.9, 0.98129), (1.0, 1.0))
        Δt = cfl*Δx
        k = 90π/L
        vp = (2/Δt*asin(min(1.0, cfl*sin(k*Δx/2))))/k
        println("  cfl = ", rpad(cfl, 4), " at kΔx = ", round(k*Δx; digits = 3),
                ": phase velocity = ", round(vp; digits = 5), "c")
        @test isapprox(vp, expected; atol = 1e-4)
    end
end

@testset "PML reflection falls with layer thickness" begin
    # The existing PML test measures the reflection coefficient at one
    # configuration, `N = 10, σ_max = 1e3`. That establishes the layer absorbs;
    # it does not establish that the σ profile is doing the absorbing, which a
    # layer that happened to be lossy at one thickness would also pass.
    #
    # Sweeping the thickness tests the ramp. Measured, σ_max = 1e3, a Gaussian
    # pulse over 800 steps:
    #
    #   N_pml    R
    #   2        8.73e-03
    #   4        2.97e-05      295x better
    #   8        2.16e-08      1374x
    #   16       3.74e-10      58x
    #   32       5.84e-12      64x
    #
    # Four orders between 2 and 8 cells is the cubic ramp working. The rate
    # falls off beyond that, which is the expected shape -- the reflection stops
    # being limited by the layer and starts being limited by the discretisation
    # of the ramp itself -- so the assertions below are strongest where the
    # physics is.
    Δx = 0.01
    cfl = 0.8
    Δt = cfl*Δx
    N = 400

    function residual(NP)
        m = FDTD1D.YeeMesh1D{Float64}(N)
        for i = 0:N
            m.ey[i+1] = exp(-((i - 200)/12)^2)
            i + 1 ≤ N && (m.hz[i+1] = exp(-((i + 0.5 - 200 - 0.5*cfl)/12)^2))
        end
        incident = maximum(abs, m.ey)
        op = FDTD1D.Yee1D(; Δx, Δt, cfl, source = NO_PULSE, pml = FDTD1D.PML(NP, 1e3, Δx, Δt))
        j = zero_current(N + 1)
        for s = 1:800
            FDTD1D.advance!(m, op, s*Δt, j)
        end
        return maximum(abs, m.ey[NP+2:N-NP])/incident
    end

    thicknesses = (2, 4, 8, 16, 32)
    R = [residual(NP) for NP in thicknesses]
    for (NP, r) in zip(thicknesses, R)
        println("  N_pml = ", lpad(NP, 3), "   R = ", r)
    end

    @test issorted(R; rev = true)             # thicker is always better
    @test R[2] < R[1]/100                     # 2 -> 4 cells: measured 295x
    @test R[3] < R[2]/100                     # 4 -> 8 cells: measured 1374x
    @test R[3] < 1e-7
    @test R[end] < 1e-10
end
