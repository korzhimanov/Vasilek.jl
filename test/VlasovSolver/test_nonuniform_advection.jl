using Vasilek: PFCNonUniform, PFC

@testset "PFCNonUniform" begin
    # ---- the limiter coefficient -------------------------------------------
    # A locally uniform triple must give exactly 2, the standard PFC value.
    @test PFCNonUniform.slope_limit(1.0) ≈ 2.0
    # and a ratio of 1:2 must tighten it well below 2
    @test PFCNonUniform.slope_limit(0.5) < 1.0

    # ---- mass conservation --------------------------------------------------
    # The scheme is written in flux form: whatever leaves cell i enters its
    # neighbour, so Σ f·Δx is conserved to round-off on a periodic grid.
    for Δx in (fill(0.1, 64), vcat(fill(0.2, 16), fill(0.1, 32), fill(0.2, 16)))
        f = [exp(-((i*0.1 - 3.2)/0.5)^2) for i = 0:63]
        m₀ = sum(f .* Δx)
        advect! = PFCNonUniform.make_advect_1D!(Δx; fₘᵢₙ = 0.0, fₘₐₓ = 1.0)
        for _ = 1:50
            advect!(f, 0.04)
        end
        @test sum(f .* Δx) ≈ m₀ rtol=1e-13
    end

    # ---- positivity and the maximum principle ------------------------------
    # This is the entire reason the scheme exists, and nothing tested it.
    Δx = vcat(fill(0.2, 16), fill(0.1, 32), fill(0.2, 16))
    f = [exp(-((i*0.1 - 3.2)/0.5)^2) for i = 0:63]
    advect! = PFCNonUniform.make_advect_1D!(Δx; fₘᵢₙ = 0.0, fₘₐₓ = 1.0)
    for _ = 1:200
        advect!(f, 0.04)
        advect!(f, -0.04)
    end
    @test minimum(f) ≥ -eps()
    @test maximum(f) ≤ 1.0 + eps()

    # ---- no allocation on the steady-state path ----------------------------
    # The previous implementation sliced f[i-1:i+1] and Δx[i-1:i+1], allocating
    # two three-element vectors per cell per step, and copied in and out.
    g = [exp(-((i*0.1 - 3.2)/0.5)^2) for i = 0:63]
    step! = PFCNonUniform.make_advect_1D!(Δx; fₘᵢₙ = 0.0, fₘₐₓ = 1.0)
    step!(g, 0.04)                       # compile
    @test (@allocated step!(g, 0.04)) == 0

    # ---- agreement with the uniform-grid solver ----------------------------
    # On a uniform grid the non-uniform scheme must reproduce PFC. Both reduce
    # to the same formulas with ξ = 2 and Δx cancelling out.
    h = [exp(-((i*0.1 - 3.2)/0.5)^2) for i = 0:63]
    nonuni = copy(h)
    step_nu! = PFCNonUniform.make_advect_1D!(fill(0.1, 64); fₘᵢₙ = 0.0, fₘₐₓ = 1.0)
    step_nu!(nonuni, 0.4*0.1)            # α is a length, so c = α/Δx = 0.4

    uni = similar(h)
    PFC.generate_solver(h, uni; fₘᵢₙ = 0.0, fₘₐₓ = 1.0)(0.4)
    println("PFCNonUniform vs PFC on a uniform grid: max|Δ| = ",
            maximum(abs, nonuni .- uni))
    @test nonuni ≈ uni rtol=1e-12
end
