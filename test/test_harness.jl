@isdefined(vlasov_poisson) || include(joinpath(@__DIR__, "verification_harness.jl"))

@testset "The harness's Poisson solve takes any grid" begin
    # Its wavenumbers used to come from a float range whose length rounded, so
    # 61 of the even Nx in 8:512 on the 8π box threw a DimensionMismatch on the
    # first solve -- every grid the suite happened to use was a lucky one.
    #
    # The harness now calls the package solver, so the reference has to be
    # independent of it: the analytic field. With de/dx = ρ, ρ = cos(k₁x) +
    # 0.3 sin(k₂x) gives e = sin(k₁x)/k₁ − 0.3 cos(k₂x)/k₂, and the centred
    # difference the harness uses scales each mode by sin(kΔx)/(kΔx). The
    # spectral form must give the field itself.
    L = 8π
    k₁, k₂ = 0.25, 0.75
    failures = 0
    worst = 0.0
    worst_spectral = 0.0
    for Nx in 8:1:512
        Δx = L/Nx
        x = collect(range(Δx; step = Δx, length = Nx))
        ρ = @. cos(k₁*x) + 0.3sin(k₂*x)
        s₁, s₂ = sin(k₁*Δx)/(k₁*Δx), sin(k₂*Δx)/(k₂*Δx)
        exact    = @. sin(k₁*x)/k₁ - 0.3cos(k₂*x)/k₂
        centered = @. s₁*sin(k₁*x)/k₁ - s₂*0.3cos(k₂*x)/k₂
        e = try
            make_poisson(x)(similar(x), ρ)
        catch err
            err isa DimensionMismatch || rethrow()
            failures += 1
            continue
        end
        worst = max(worst, maximum(abs, e .- centered))
        p = PoissonFourier1D.PoissonFFT1D(Nx, Δx; derivative = :spectral)
        worst_spectral = max(worst_spectral,
                             maximum(abs, PoissonFourier1D.solve!(similar(x), ρ, p) .- exact))
    end
    println("  make_poisson over Nx = 8:512: ", failures, " failures, max|Δe| vs analytic = ",
            worst, " (centred), ", worst_spectral, " (spectral)")
    @test failures == 0
    @test worst < 1e-12
    @test worst_spectral < 1e-12
end
