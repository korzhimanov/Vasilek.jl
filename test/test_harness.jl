@isdefined(vlasov_poisson) || include(joinpath(@__DIR__, "verification_harness.jl"))

@testset "The harness's Poisson solve takes any grid" begin
    # Its wavenumbers used to come from a float range whose length rounded, so
    # 61 of the even Nx in 8:512 on the 8π box threw a DimensionMismatch on the
    # first solve -- every grid the suite happened to use was a lucky one. The
    # answer must also be the package solver's, which uses the same spectrum.
    L = 8π
    failures = 0
    worst = 0.0
    for Nx in 8:1:512
        Δx = L/Nx
        x = collect(range(Δx; step = Δx, length = Nx))
        ρ = @. cos(0.25x) + 0.3sin(0.75x)
        e = try
            make_poisson(x)(similar(x), ρ)
        catch err
            err isa DimensionMismatch || rethrow()
            failures += 1
            continue
        end
        ref = PoissonFourier1D.solve!(similar(x), ρ, PoissonFourier1D.PoissonFFT1D(Nx, x[2] - x[1]))
        worst = max(worst, maximum(abs, e .- ref))
    end
    println("  make_poisson over Nx = 8:512: ", failures, " failures, max|Δe| vs package = ", worst)
    @test failures == 0
    @test worst < 1e-12
end
