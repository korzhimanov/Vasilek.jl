# Extended verification. About 40 seconds now that PFCNonUniform's per-cell
# allocations are gone -- cheap enough to run far more often than nightly --
# but still gated so the default test run stays instant:
#
#     VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'
#
# These promote the claims the verification notebooks make in prose into
# assertions. Until they run, those claims are a 2021 HTML file.

include(joinpath(@__DIR__, "verification_harness.jl"))

@testset "Extended verification" begin
    if get(ENV, "VASILEK_EXTENDED", "0") != "1"
        @info "extended verification skipped; set VASILEK_EXTENDED=1 to run it"
        @test true
    else
        @testset "Landau damping, k = 0.5" begin
            x = collect(π/8:π/8:8π)
            v = collect(-4:0.1:4)
            f₀ = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(0.5*x)))'
            t = collect(0.0:0.1:70.0)
            ε_e, _ = vlasov_poisson(x, v, f₀, t)

            # The electric energy decays as exp(-2γt). Fit over the clean
            # stretch before recurrence sets in.
            i0, i1 = 60, 300
            A = hcat(ones(i1 - i0 + 1), t[i0:i1])
            fitted = -(A \ log.(ε_e[i0:i1]))[2] / 2

            # Canosa's tabulated root for kλ_D = 0.5. Note this is *not* the
            # π/(8√2)u³exp(-u²/2) asymptotic the notebook quotes, which gives
            # 0.1145 and is not accurate at this k -- measuring against it
            # suggested a 31% error where there is none.
            analytic = 0.1533
            println("  Landau γ: fitted ", round(fitted; digits = 4),
                    " vs analytic ", analytic,
                    "  (", round(100*abs(fitted - analytic)/analytic; digits = 1), "%)")
            @test isapprox(fitted, analytic; rtol = 0.05)
        end

        @testset "Plasma oscillations, energy conservation" begin
            x = collect(1.0:1.0:100.0)
            t = collect(0.0:0.1:3000.0)
            f(v) = 1/sqrt(2π)*(@. exp(-0.5*v^2)) * (@. (1.0 + 0.01*cos(2π*x/100)))'

            # The uniform-grid claim the notebook makes: below 0.5% at t = 3000.
            uniform = collect(-4:0.1:4)
            _, ε = vlasov_poisson(x, uniform, f(uniform), t)
            drift = (ε[end-1] - ε[1])/ε[1]
            println("  uniform grid Δε/ε = ", drift)
            @test abs(drift) < 0.005

            # The non-uniform grid was reported at "about 12%" with the limiter
            # computed globally. Per-cell-triple ξ brings it to 4.85%.
            nonuniform = vcat(collect(-4:0.2:-1.2), collect(-1:0.1:1), collect(1.2:0.2:4))
            _, ε = vlasov_poisson(x, nonuniform, f(nonuniform), t)
            drift = (ε[end-1] - ε[1])/ε[1]
            println("  non-uniform grid Δε/ε = ", drift)
            @test abs(drift) < 0.06
        end
    end
end
