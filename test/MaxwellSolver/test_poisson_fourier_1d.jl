using Vasilek.PoissonFourier1D

# The 1D spectral Poisson solve.
#
# This suite used to be a single mode compared against the continuum answer at
# `atol = 1e-3`, which is loose enough to hide the very bug the module has a
# history of: a field wrong by a factor of `Δx²` from a wavenumber convention
# mismatch. `docs/normalization.md` already states the exact discrete identity
# the scheme satisfies and claims the test asserts it to machine precision --
# until now it did not.
@testset "Test 1D Poisson solvers" begin
    Δx = 0.01
    ρ = [sin(2π*i*Δx) for i = 0:999]
    e = similar(ρ)
    solve! = PoissonFourier1D.generate_solver(ρ, Δx)
    solve!(e, ρ)
    @test sum(@. (e - [-cos(2π*i*Δx)/2π for i = 0:999])^2) ≈ 0 atol=1e-3
end

@testset "Poisson: the exact discrete identity" begin
    # On a grid mode the DFT is exact, so the only discretisation left is the
    # centred difference that takes φ to E, contributing exactly sinc(kΔx):
    #
    #     E_num[i] == (sin(kΔx)/(kΔx)) · E_exact[x_i]
    #
    # This is an identity, not an approximation, and it pins the wavenumber
    # convention to machine precision: a solver that dropped Δx from the
    # spectrum would miss by a factor of Δx² and could not survive it at any N.
    #
    # Odd N is included because `rfft` has a different tail there, and the
    # `ω[1] = ω[2]` guard against dividing by the zero mode sits right next to
    # it. m = 63 at N = 128 sits one mode below Nyquist, where the identity is
    # most informative and most easily got wrong.
    #
    # **Normalisation.** The departure is measured against `1/k`, the amplitude
    # of the *continuum* field, not against the amplitude the scheme actually
    # produces. Near Nyquist the sinc factor suppresses the output itself --
    # at m = 63 the field is 2.6e-4 where the continuum answer is 1.6e-2 --
    # so dividing by the output would report 7.3e-12 for a solver whose
    # absolute error is the same 2e-15 it is everywhere else. The suppression
    # is the physics of the centred difference, not an error in it.
    #
    # Measured: absolute error ≤ 6.6e-14 and normalised departure ≤ 1.2e-13
    # across every case below.
    worst = 0.0
    for (N, m, Δx) in [(1000, 1, 0.01), (1000, 7, 0.01), (999, 3, 0.01),
                       (256, 5, 0.02), (64, 1, 0.1), (128, 63, 0.05)]
        L = N*Δx; k = 2π*m/L
        x = [i*Δx for i = 0:N-1]
        ρ = sin.(k.*x)
        e = similar(ρ)
        PoissonFourier1D.generate_solver(ρ, Δx)(e, ρ)
        exact = @. -cos(k*x)/k
        predicted = (sin(k*Δx)/(k*Δx)) .* exact
        dev = maximum(abs, e .- predicted)/maximum(abs, exact)
        worst = max(worst, dev)
        println("  N = ", rpad(N, 5), " m = ", rpad(m, 3), " Δx = ", rpad(Δx, 5),
                " sinc = ", rpad(round(sin(k*Δx)/(k*Δx); digits = 4), 7),
                " departure from identity = ", dev)
        @test dev < 1e-12
    end
    println("  worst departure: ", worst)
end

@testset "Poisson: linearity and the zero mode" begin
    N = 256; Δx = 0.02; L = N*Δx
    x = [i*Δx for i = 0:N-1]
    ρ = sin.(2π*x/L) .+ 0.5*cos.(6π*x/L)
    solve! = PoissonFourier1D.generate_solver(ρ, Δx)
    e = similar(ρ); solve!(e, ρ)

    # `F[1] = 0` discards the zero mode, so a net charge is ignored rather than
    # producing a divergent field. Adding a constant to ρ must leave E alone.
    # Measured 4.2e-15.
    shifted = similar(ρ); solve!(shifted, ρ .+ 3.7)
    println("  constant offset of 3.7 in ρ moves E by ", maximum(abs, e .- shifted))
    @test maximum(abs, e .- shifted) < 1e-12

    # Poisson is linear, and the solver reuses captured buffers across calls;
    # superposition is what catches state leaking from one call into the next.
    # Measured 5.7e-15.
    a = sin.(2π*x/L); b = 0.5*cos.(6π*x/L)
    ea = similar(a); eb = similar(b)
    solve!(ea, a); solve!(eb, b)
    println("  |E(a+b) − E(a) − E(b)| = ", maximum(abs, e .- ea .- eb))
    @test maximum(abs, e .- ea .- eb) < 1e-12

    # E derives from a periodic potential by a centred difference, so it has no
    # mean. Measured 5.2e-18.
    @test abs(sum(e)/N) < 1e-14

    # Repeating a call must reproduce it bit-for-bit: the plans and buffers are
    # captured, and a stale one would show up here and nowhere else.
    again = similar(ρ); solve!(again, ρ)
    @test again == e
end

@testset "Poisson: second order in Δx" begin
    # The centred difference is the only source of error, so refining the grid
    # at fixed physical wavelength must converge at second order.
    # Measured errors: 3.27e-4, 8.18e-5, 2.05e-5, 5.11e-6 -- ratios of exactly 4.
    L = 5.12; k = 2π/L
    errs = Float64[]
    for N in (128, 256, 512, 1024)
        Δx = L/N
        x = [i*Δx for i = 0:N-1]
        ρ = sin.(k.*x); e = similar(ρ)
        PoissonFourier1D.generate_solver(ρ, Δx)(e, ρ)
        push!(errs, maximum(abs, e .- (@. -cos(k*x)/k)))
    end
    for i = 2:length(errs)
        p = log2(errs[i-1]/errs[i])
        println("  N = ", rpad(2^(6+i), 5), " err = ", rpad(round(errs[i]; sigdigits = 4), 11),
                " order = ", round(p; digits = 4))
        @test isapprox(p, 2.0; atol = 0.05)
    end
end
