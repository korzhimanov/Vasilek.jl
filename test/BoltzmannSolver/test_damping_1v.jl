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

"Discrete number density, mean velocity and temperature of `f` on grid `v`."
function moments(v, f)
    n = integrate(v, f)
    u = integrate(v, v.*f)/n
    return n, u, integrate(v, (v .- u).^2 .*f)/n
end

@testset "BGK conserves its moments" begin
    # Relaxation towards the *local* Maxwellian is defined by leaving n, u and T
    # alone -- it is the whole content of the operator, and nothing checked it.
    # The existing tests all start from data symmetric in v, so `u` was zero
    # throughout and the mean-velocity computation was never exercised at all.
    #
    # Measured over 100 steps on v ∈ [-10, 10], Δv = 0.05: machine precision for
    # the skewed and drifting cases, and 2.5e-10 for the bi-Maxwellian, whose
    # relaxed state is wide enough that the grid truncates its tails.
    v = collect(-10:0.05:10)
    cases = [("skewed",                  @. exp(-v^2)*(1.0 + 0.3*v^3)),
             ("drifting non-Maxwellian", @. exp(-(v-0.8)^2)*(1.0 + 0.4*(v-0.8)^2)),
             ("flat-top",                @. exp(-(v/2.2)^8)),
             ("bi-Maxwellian",           @. exp(-(v-1.0)^2) + 0.6*exp(-(v+1.5)^2/0.5))]
    for (name, f₀) in cases
        n₀, u₀, T₀ = moments(v, f₀)
        f = relax(BGK(1e-1), v, 0.1, f₀, 100)
        n₁, u₁, T₁ = moments(v, f)
        println("  ", rpad(name, 24), " Δn/n = ", abs(n₁-n₀)/abs(n₀),
                "  Δu = ", abs(u₁-u₀), "  ΔT/T = ", abs(T₁-T₀)/abs(T₀))
        @test abs(n₁ - n₀)/abs(n₀) < 1e-9
        @test abs(u₁ - u₀) < 1e-9
        @test abs(T₁ - T₀)/abs(T₀) < 1e-8
        @test minimum(f) ≥ 0.0
    end
end

@testset "Any Maxwellian is a fixed point" begin
    # Not just the symmetric unit-temperature one the suite happened to use.
    # Measured after 100 steps on v ∈ [-10, 10]: 5.6e-17 at (u, T) = (0, 0.5),
    # 1.7e-16 at (0.7, 0.5), 3.0e-15 at (-1.3, 1.0), 2.8e-17 at (1.5, 0.3).
    v = collect(-10:0.05:10)
    for (u₀, T₀) in [(0.0, 0.5), (0.7, 0.5), (-1.3, 1.0), (1.5, 0.3)]
        f₀ = @. 1/sqrt(2π*T₀)*exp(-(v - u₀)^2/(2T₀))
        f = relax(BGK(1e-2), v, 0.1, f₀, 100)
        println("  u = ", rpad(u₀, 5), " T = ", rpad(T₀, 4), " max|f - f₀| = ",
                maximum(abs, f .- f₀))
        @test maximum(abs, f .- f₀) < 1e-13
    end

    # The residue that remains is the velocity window, not the operator. A
    # Maxwellian at T = 2 has σ = 1.41, so ±8 is under six σ and the trapezoid
    # loses the tail; widening the grid recovers four orders of magnitude
    # (1.6e-5 → 3.7e-9). Quantified rather than left as a footnote, because it
    # sets the window every caller of this operator needs.
    T₀, u₀ = 2.0, 0.4
    dev(hi) = let v = collect(-hi:0.05:hi)
        f₀ = @. 1/sqrt(2π*T₀)*exp(-(v - u₀)^2/(2T₀))
        maximum(abs, relax(BGK(1e-2), v, 0.1, f₀, 100) .- f₀)
    end
    narrow, wide = dev(8.0), dev(10.0)
    println("  T = 2 truncation: window ±8 gives ", narrow, ", ±10 gives ", wide)
    @test wide < narrow/100
end

@testset "BGK limits and the exact update" begin
    # `dest = src·e + (1−e)·M` with `e = exp(-Δt/τ)`, so both limits and the
    # interpolation between them are algebraic identities, and hold exactly.
    v = collect(-10:0.05:10)
    f₀ = @. exp(-(v - 0.5)^2)*(1.0 + 0.3*v^2)
    dest = similar(f₀)
    n, u, T = moments(v, f₀)
    M = @. n/sqrt(2π*T)*exp(-(v - u)^2/(2T))

    # Δt ≪ τ: nothing happens. Measured 3.5e-12 at Δt/τ = 1e-10.
    collide!(dest, f₀, BGK(1.0), v, 1e-10, collision_workspace(BGK(1.0), length(v)))
    println("  Δt/τ = 1e-10: max|dest − src| = ", maximum(abs, dest .- f₀))
    @test maximum(abs, dest .- f₀) < 1e-10

    # Δt ≫ τ: the local Maxwellian, bit-for-bit -- `e` underflows to zero and
    # the update collapses to `M` exactly.
    collide!(dest, f₀, BGK(1e-8), v, 1.0, collision_workspace(BGK(1e-8), length(v)))
    @test dest == M
    @test minimum(dest) ≥ 0.0

    # and in between, the stated formula, also bit-for-bit
    τ, Δt = 0.7, 0.3
    collide!(dest, f₀, BGK(τ), v, Δt, collision_workspace(BGK(τ), length(v)))
    e = exp(-Δt/τ)
    @test dest == @. f₀*e + (1.0 - e)*M
end
