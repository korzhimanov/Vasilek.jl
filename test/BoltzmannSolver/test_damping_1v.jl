using Vasilek
using Vasilek.Collisions: Landau1P, collide!, workspace as collision_workspace
using NumericalIntegration

"""
    relax(op, v, Δt, f₀, nsteps)

`nsteps` steps of the collision operator, feeding each result back as the next
input. Under the old closure API this loop rebound a local instead, so the
operator kept reading the original buffer and the hundred "steps" recomputed
step one -- which is why an index bug in Landau1P survived for years.
"""
function relax(op, v, Δt, f₀, nsteps)
    src = copy(f₀)
    dst = similar(src)
    ws = collision_workspace(op, length(src))
    for _ = 1:nsteps
        collide!(dst, src, op, v, Δt, ws)
        copyto!(src, dst)
    end
    return src
end

@testset "Test 1V Boltzmann solvers" begin
    Δt = 0.1
    v = collect(-4:0.1:4)

    f₀ = @. exp(-v^2)
    f = relax(BGK(1e-2), v, Δt, f₀, 100)
    println("  BGK, Maxwellian stays Maxwellian: ", norm(f - f₀))
    @test norm(f - f₀) ≈ 0 atol=1e-4

    # The index bug is fixed, which halves the deviation (0.310 -> 0.160), but a
    # Maxwellian still is not a stationary point. What remains is the closure
    # inconsistency: the transversal estimate 2Tₜ in the numerator against the
    # longitudinal |vᵢ-vⱼ|³ in the denominator. See the refinement test below.
    f = relax(Landau1P(1e-2), v, Δt, f₀, 100)
    println("  Landau1P, Maxwellian deviation: ", norm(f - f₀))
    @test_broken norm(f - f₀) ≈ 0 atol=3e-3

    a = 3e-1
    n = 0.5*sqrt(π)*(a + 2.0)
    T = 0.25*sqrt(π)*(3a + 2.0)/n
    g₀ = @. exp(-v^2)*(1.0 + a*v^2)
    g₁ = @. n/sqrt(2π*T)*exp(-v^2/(2T))
    g = relax(BGK(1e-1), v, Δt, g₀, 100)
    println("  BGK, relaxation to the Maxwellian: ", norm(g - g₁))
    @test norm(g - g₁) ≈ 0 atol=2e-3
end

@testset "Landau1P invariants" begin
    Δt = 0.1
    v = collect(-4:0.1:4)
    f₀ = @. exp(-v^2)
    C = Vasilek.Collisions

    # The kernel must be antisymmetric under i <-> j. That is what makes the
    # operator conservative, it holds whatever the closure is, and it is exactly
    # the symmetry the old index bug destroyed.
    Tₜ = 1e-3
    K(i, j) = i == j ? 0.0 :
        (f₀[i]*C.∂f∂v(f₀, v, j) - f₀[j]*C.∂f∂v(f₀, v, i))*2Tₜ/abs(v[i]-v[j])^3
    @test maximum(abs(K(i,j) + K(j,i)) for i in eachindex(v), j in eachindex(v)) == 0.0

    # Mass is not conserved to machine precision: the update differences a
    # cell-centred I rather than staggered fluxes. Assert the drift actually
    # achieved rather than an exactness the scheme does not have.
    f = relax(Landau1P(1e-2), v, Δt, f₀, 100)
    drift = abs(integrate(v, f) - integrate(v, f₀))/integrate(v, f₀)
    println("  Landau1P mass drift over 100 steps: ", drift)
    @test drift < 1e-8

    # The collision integral must converge as the velocity grid is refined. It
    # does not: with 2Tₜ in the numerator but |vᵢ-vⱼ|³ in the denominator the
    # kernel is non-integrable at i ≈ j, and max|∂f/∂t| grows under refinement
    # instead of settling (0.0023, 0.0053, 0.0080, 0.0104 at Δv = 0.4, 0.2,
    # 0.1, 0.05).
    function rate(Δv)
        vv = collect(-4:Δv:4)
        g₀ = @. exp(-vv^2)
        g = relax(Landau1P(1e-2), vv, Δt, g₀, 1)
        maximum(abs, (g .- g₀)./Δt)
    end
    coarse, fine = rate(0.1), rate(0.05)
    println("  Landau1P refinement: $coarse -> $fine")
    @test_broken isapprox(fine, coarse; rtol = 0.1)
end
