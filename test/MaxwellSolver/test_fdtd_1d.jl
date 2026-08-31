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
        @test length(m.ey) == length(m.ez) == 8
        @test length(m.hy) == length(m.hz) == 7
        @test m.N == 7
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
        advance! = FDTD1D.make_advance_fields(m, cfl, NO_PULSE, Δt, Δx, 0, no_pml(Δx, Δt))
        j = zero_current(N + 1)

        staggered = Float64[]; equivalent = Float64[]; naive = Float64[]
        for s = 1:2000
            e_pre = copy(m.ey); h_pre = copy(m.hz)
            advance!(s*Δt, j)
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
        advance! = FDTD1D.make_advance_fields(m, cfl, NO_PULSE, Δt, Δx, 0, pml)
        j = zero_current(N + 1)
        for s = 1:600
            advance!(s*Δt, j)
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
        advance! = FDTD1D.make_advance_fields(m, cfl, NO_PULSE, Δt, Δx, 0, pml)
        j = zero_current(N + 1)
        for s = 1:50
            advance!(s*Δt, j)
        end
        return m
    end
    bare = run(no_pml(Δx, Δt))
    faint = run(FDTD1D.PML(10, 1.25e-4, Δx, Δt))
    dev = maximum(abs, faint.ey[15:N-14] .- bare.ey[15:N-14])
    println("  σ_max = 1.25e-4 vs no layer: max|Δ| in the interior = ", dev)
    @test dev < 1e-12
end

@testset "Where the current is applied, and which nodes are frozen" begin
    # Two asymmetries in the field update, pinned as they stand rather than
    # changed -- both are boundary conventions, and which one is intended is the
    # author's call.
    #
    #   * `update_ey!` adds the current over `1:Nx`, but `ey` has `Nx+1` entries
    #     and every caller in this repository passes a current of length `Nx+1`.
    #     The last node's current is silently dropped.
    #   * Neither `ey[1]` nor `ey[end]` is ever touched by a curl loop, at any
    #     PML setting: the interior loop runs `pml.N+2 : Nx-pml.N` and the two
    #     layer loops stop short of both. They are frozen at their initial
    #     values, which is a PEC boundary by omission -- and `ey[end]` is read
    #     by `update_hz!`, so it is a real boundary condition rather than dead
    #     storage.
    #
    # Together: `ey[1]` receives current but no curl, `ey[end]` receives
    # neither. If the intent is PEC at both ends, the current at node 1 is the
    # inconsistent one.
    Δx = 0.01; cfl = 0.8; Δt = cfl*Δx; N = 20
    for NP in (0, 5)
        m = FDTD1D.YeeMesh1D{Float64}(N)
        advance! = FDTD1D.make_advance_fields(m, cfl, NO_PULSE, Δt, Δx, 0,
                                              FDTD1D.PML(NP, 1e3, Δx, Δt))
        advance!(0.0, (y = fill(1.0, N+1), z = fill(1.0, N+1)))
        @test length(m.ey) == N + 1
        @test all(m.ey[1:N] .== 1.0)        # current applied here
        @test m.ey[N+1] == 0.0              # and dropped here
    end

    for NP in (0, 10)
        m = FDTD1D.YeeMesh1D{Float64}(60)
        m.ey .= 7.0; m.hz .= 3.0
        advance! = FDTD1D.make_advance_fields(m, cfl, NO_PULSE, Δt, Δx, 0,
                                              FDTD1D.PML(NP, 1e3, Δx, Δt))
        j = zero_current(61)
        for s = 1:5
            advance!(s*Δt, j)
        end
        @test m.ey[1] == 7.0                # frozen: no curl loop reaches either end
        @test m.ey[end] == 7.0
    end
end
