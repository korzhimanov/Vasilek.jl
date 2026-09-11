# The two-stream instability: growth rate against the cold dispersion relation.
#
#     julia --project=verification verification/two-stream.jl
#
# Writes two-stream-*.png beside this script.
#
# Every other verification here is a damped or neutral mode. This is the first
# *unstable* one, and it catches a class of error damping cannot: a sign flip
# in the field push turns damping into growth and growth into damping, so a
# suite made only of damped cases is half-blind to it.
#
# `two_stream`, `γ_cold` and friends mirror the definitions in
# `test/test_verification.jl` rather than living in `test/verification_harness.jl`
# -- the same choice `landau_case` makes there, and for the same reason: they are
# the *setup* for one problem, not a piece of the shared driver. `vlasov_poisson`,
# `cell_widths` and `growth_rate` are the parts that are shared, and come from the
# harness below so that this script and the test measure the same solver.

using Plots

include(joinpath(@__DIR__, "..", "test", "verification_harness.jl"))

here = @__DIR__

"""
    γ_cold(a)

Growth rate of the cold two-stream instability at `a = kv₀`, for two beams of
density 1/2 at `±v₀` with `ω_p = 1`. The electrostatic dispersion relation

    1 = ½/(ω − kv₀)² + ½/(ω + kv₀)²

is, with `u = ω²`, a quadratic `u² − (2a² + 1)u + (a⁴ − a²) = 0`, and `u₋ < 0`
exactly when `a < 1` -- then `γ = √(−u₋)`. Closed form, so this needs no
tabulated constant of the kind the Landau cases carry.
"""
two_stream_u(a) = ((2a^2 + 1) - sqrt(8a^2 + 1))/2
γ_cold(a) = sqrt(max(0.0, -two_stream_u(a)))

"The cold dispersion relation itself, for checking `γ_cold` against."
two_stream_residual(ω, a) = 0.5/(ω - a)^2 + 0.5/(ω + a)^2 - 1

"""
    two_stream(a; v₀, vt, Δv, vmax, Δt, tmax)

Two counter-streaming warm beams at `±v₀`, perturbed by 0.1% in the `k = a/v₀`
mode, returning `(t, ε_e)` over one wavelength.
"""
function two_stream(a; v₀ = 3.0, vt = 0.3, Δv = 0.05, vmax = 6.0,
                       Δt = 0.05, tmax = 24.0)
    k = a/v₀
    L = 2π/k
    Nx = round(Int, L/0.49)
    Δx = L/Nx
    x = collect(Δx:Δx:L)
    v = collect(-vmax:Δv:vmax)
    t = collect(0.0:Δt:tmax)
    beams = @. 0.5/sqrt(2π*vt^2)*(exp(-(v - v₀)^2/(2vt^2)) +
                                  exp(-(v + v₀)^2/(2vt^2)))
    f₀ = beams * (@. (1.0 + 1e-3*cos(k*x)))'
    # The harness defaults to `fmax = 1.0`, which a beam this narrow exceeds:
    # the peak is 0.5/(√(2π)·vt) = 0.665 at vt = 0.3, and passes 1.0 below 0.2.
    r = vlasov_poisson(x, v, f₀, t;
            scheme_x = PFCNonUniform(cell_widths(x); fmin = 0.0, fmax = 3.0),
            scheme_v = PFCNonUniform(cell_widths(v); fmin = 0.0, fmax = 3.0))
    return t[1:end-1], r.ε_e[1:end-1]
end

# Sanity check before trusting the closed form for the rest of the script.
worst = maximum(a -> abs(two_stream_residual(im*γ_cold(a), a)), (0.2, 0.4, 0.6, 0.8, 0.95))
println("closed form vs dispersion relation: worst |residual| = ", worst)

# ---- growth curves for the three wavenumbers the test asserts, with the
# fitted exponential overlaid on the window `growth_rate` actually used.
growth = plot(yscale = :log10, xlabel = "ωₚt", ylabel = "εₑ",
              legend = :bottomright,
              title = "Two-stream instability: field-energy growth",
              size = (700, 460))
measured_a = (0.4, 0.6, 0.8)
measured_γ = Float64[]
for (a, color) in zip(measured_a, (:steelblue, :crimson, :seagreen))
    t, ε_e = two_stream(a)
    γ, t0, t1 = growth_rate(t, ε_e; lo = 100*ε_e[1], hi = 5.0)
    push!(measured_γ, γ)
    plot!(growth, t, ε_e; label = "a = $a  (γ = $(round(γ; digits = 3)))",
          color = color, linewidth = 1.8)
    i0 = argmin(abs.(t .- t0))
    tt = range(t0, t1; length = 50)
    plot!(growth, tt, ε_e[i0] .* exp.(γ .* (tt .- t0));
          linestyle = :dash, color = color, label = "")
end
savefig(growth, joinpath(here, "two-stream-growth.png"))

# ---- the growth rate against the closed form over the whole branch: γ(a) is
# non-monotone, peaking at a = √(3/8) and reaching zero at the stability
# boundary a = 1 -- reproducing that shape is a statement about the dispersion
# relation, not about one point on it.
avals = 0.001:0.002:1.25
dispersion = plot(avals, γ_cold.(avals); label = "γ_cold(a)  (closed form)",
                   xlabel = "a = kv₀", ylabel = "γ",
                   title = "Two-stream growth rate vs wavenumber",
                   linewidth = 2.2, color = :steelblue, size = (700, 460))
scatter!(dispersion, collect(measured_a), measured_γ;
         label = "measured (vt = 0.3)", markersize = 6, color = :crimson)
vline!(dispersion, [sqrt(3/8)]; linestyle = :dot, color = :gray,
       label = "a = √(3/8) (peak)")
vline!(dispersion, [1.0]; linestyle = :dashdot, color = :black,
       label = "a = 1 (stability boundary)")
savefig(dispersion, joinpath(here, "two-stream-dispersion.png"))

println("a, γ measured, γ cold, error")
for (a, γ) in zip(measured_a, measured_γ)
    println("  ", a, "  ", round(γ; digits = 5), "  ", round(γ_cold(a); digits = 5),
            "  (", round(100*(γ - γ_cold(a))/γ_cold(a); digits = 2), "%)")
end
println("wrote two-stream-{growth,dispersion}.png to ", here)
