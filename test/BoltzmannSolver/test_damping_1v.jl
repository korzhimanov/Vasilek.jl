using Vasilek: BGK, Landau1P

"""
    test_1v_damping_step(mod, Δt, f₀, f₁, v, exp_norm_dev, args...; kwargs...)

Relax `f₀` for 100 steps of `mod` and compare against `f₁`.

`Landau1P` builds a closure over the *arrays* it was handed, so advancing the
state means writing into that same array — rebinding the local name leaves the
closure reading the original buffer and silently repeats step one forever.
`BGK` instead mutates its single buffer in place.
"""
function test_1v_damping_step(mod, Δt, f₀, f₁, v, exp_norm_dev, args...; kwargs...)
    name = nameof(mod)

    if mod === BGK
        f = copy(f₀)
        damp! = mod.generate_solver(f, v, Δt, args...; kwargs...)
        for _ = 1:100
            damp!()
        end
    else
        src = copy(f₀)
        f = similar(f₀)
        damp! = mod.generate_solver(src, f, v, Δt, args...; kwargs...)
        for _ = 1:100
            damp!()
            copyto!(src, f)   # feed the result back through the captured buffer
        end
    end

    dev = norm(f - f₁)
    println("$name $args $(values(values(kwargs))): $dev")
    return f, dev
end

@testset "Test 1V Boltzmann solvers" begin
    Δt = 0.1
    v = collect(-4:0.1:4)

    f₀ = @. exp(-v^2)
    f₁ = @. exp(-v^2)

    _, dev = test_1v_damping_step(BGK, Δt, f₀, f₁, v, 1e-4, 1e-2)
    @test dev ≈ 0 atol=1e-4

    # The index bug is fixed, which halves the deviation (0.310 -> 0.160), but a
    # Maxwellian still is not a stationary point. What remains is the closure
    # inconsistency: the transversal estimate 2Tₜ in the numerator against the
    # longitudinal |vᵢ-vⱼ|³ in the denominator. See the refinement test below.
    _, dev = test_1v_damping_step(Landau1P, Δt, f₀, f₁, v, 3e-3, 1e-2)
    @test_broken dev ≈ 0 atol=3e-3

    a = 3e-1
    n = 0.5*sqrt(π)*(a+2.)
    T = 0.25*sqrt(π)*(3a+2.)/n
    f₀ = @. exp(-v^2)*(1.0 + a*v^2)
    f₁ = @. n/sqrt(2π*T)*exp(-v^2/(2T))

    _, dev = test_1v_damping_step(BGK, Δt, f₀, f₁, v, 2e-3, 1e-1)
    @test dev ≈ 0 atol=2e-3
end

@testset "Landau1P invariants" begin
    using NumericalIntegration

    Δt = 0.1
    v = collect(-4:0.1:4)
    f₀ = @. exp(-v^2)

    # The kernel must be antisymmetric under i <-> j. That is what makes the
    # operator conservative, and it holds whatever the closure is. It is also
    # exactly the symmetry the old index bug destroyed: branching the inner
    # derivative on the outer index made Δfⱼ ≡ Δfᵢ.
    Tₜ = 1e-3
    K(i, j) = i == j ? 0.0 :
        (f₀[i]*Landau1P.∂f∂v(f₀, v, j) - f₀[j]*Landau1P.∂f∂v(f₀, v, i)) *
        2Tₜ/abs(v[i]-v[j])^3
    @test maximum(abs(K(i,j) + K(j,i)) for i in eachindex(v), j in eachindex(v)) == 0.0

    # Mass is not conserved to machine precision: the update differences a
    # cell-centred I rather than staggered fluxes. Assert the drift actually
    # achieved (2e-10 over 100 steps) so a regression is visible, rather than
    # asserting an exactness the scheme does not have.
    src = copy(f₀)
    f = similar(f₀)
    step! = Landau1P.generate_solver(src, f, v, Δt, 1e-2)
    m₀ = integrate(v, f₀)
    for _ = 1:100
        step!()
        copyto!(src, f)
    end
    drift = abs(integrate(v, f) - m₀)/m₀
    println("Landau1P mass drift over 100 steps: $drift")
    @test drift < 1e-8

    # The collision integral must converge as the velocity grid is refined.
    # It does not: with the transversal estimate 2Tₜ in the numerator but the
    # longitudinal |vᵢ-vⱼ|³ in the denominator, the kernel is non-integrable
    # at i ≈ j, and max|∂f/∂t| grows under refinement instead of settling
    # (measured 0.0023, 0.0053, 0.0080, 0.0104 at Δv = 0.4, 0.2, 0.1, 0.05).
    rate(Δv) = begin
        vv = collect(-4:Δv:4)
        g₀ = @. exp(-vv^2)
        s = copy(g₀); g = similar(g₀)
        Landau1P.generate_solver(s, g, vv, Δt, 1e-2)()
        maximum(abs, (g .- s)./Δt)
    end
    coarse, fine = rate(0.1), rate(0.05)
    println("Landau1P refinement: $coarse -> $fine")
    @test_broken isapprox(fine, coarse; rtol = 0.1)
end
