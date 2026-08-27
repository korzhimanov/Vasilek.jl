using Vasilek

include(joinpath(@__DIR__, "..", "scheme_cases.jl"))

@testset "PFCNonUniform" begin
    # ---- the limiter coefficient -------------------------------------------
    A = Vasilek.Advection
    @test A.slope_limit(1.0) ≈ 2.0        # locally uniform triple
    @test A.slope_limit(0.5) < 1.0        # a 1:2 ratio tightens it well below 2

    uniform = fill(0.1, 64)
    refined = vcat(fill(0.2, 16), fill(0.1, 32), fill(0.2, 16))
    f₀ = [exp(-((i*0.1 - 3.2)/0.5)^2) for i = 0:63]

    # ---- mass conservation --------------------------------------------------
    # Flux form: whatever leaves cell i enters its neighbour, so Σ f·Δx is
    # conserved to round-off on a periodic grid.
    for Δx in (uniform, refined)
        scheme = PFCNonUniform(Δx; fmin = 0.0, fmax = 1.0)
        m₀ = sum(f₀ .* Δx)
        f = march!(f₀, scheme, 0.04, 50)
        @test sum(f .* Δx) ≈ m₀ rtol=1e-13
    end

    # ---- positivity and the maximum principle ------------------------------
    # The entire reason the scheme exists, and nothing tested it before.
    scheme = PFCNonUniform(refined; fmin = 0.0, fmax = 1.0)
    ws = workspace(scheme, 64)
    src = copy(f₀); dst = similar(src)
    for _ = 1:200
        advect!(dst, src, scheme, 0.04, ws); copyto!(src, dst)
        advect!(dst, src, scheme, -0.04, ws); copyto!(src, dst)
    end
    @test minimum(src) ≥ -eps()
    @test maximum(src) ≤ 1.0 + eps()

    # ---- agreement with the uniform-grid solver ----------------------------
    # On a uniform grid the non-uniform scheme must reproduce PFC: both reduce
    # to the same formulas with ξ = 2 and Δx cancelling.
    nonuni = march!(f₀, PFCNonUniform(uniform; fmin = 0.0, fmax = 1.0), 0.4*0.1, 1)
    uni = march!(f₀, PFC(fmin = 0.0, fmax = 1.0), 0.4, 1)
    println("  PFCNonUniform vs PFC on a uniform grid: max|Δ| = ", maximum(abs, nonuni .- uni))
    @test nonuni ≈ uni rtol=1e-12

    # ---- the scheme value carries the grid, not the data -------------------
    # ξ is precomputed from Δx at construction, so one value can be applied to
    # any data on that grid, and to several tasks at once.
    other = [0.5 + 0.4*sin(2π*i/64) for i = 0:63]
    shared = PFCNonUniform(refined; fmin = 0.0, fmax = 1.0)
    @test march!(f₀, shared, 0.04, 3) == march!(f₀, PFCNonUniform(refined; fmin = 0.0, fmax = 1.0), 0.04, 3)
    @test march!(other, shared, 0.04, 3) == march!(other, PFCNonUniform(refined; fmin = 0.0, fmax = 1.0), 0.04, 3)
end
