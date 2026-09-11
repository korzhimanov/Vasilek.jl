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
# `two_stream`, `γ_cold` and the rest come from `test/verification_harness.jl`,
# not from a copy here. They were copied once, and the copy carried `tmax = 24.0`
# without the paragraph explaining that it sits in a window of about [22.9,
# 24.2] -- bounded below by the slowest fit completing and above by the fastest
# run diverging. A bare number that cannot move, with nothing saying so, is
# worse than the duplication it came from. Shared for the reason `wakefield` is
# shared: the script and the test that asserts its claims run the same setup.

using Plots

include(joinpath(@__DIR__, "..", "test", "verification_harness.jl"))

here = @__DIR__

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
    # `exp(2γΔt)`, not `exp(γΔt)`: `growth_rate` returns the rate of the field
    # *amplitude*, defined by `ε_e ∝ exp(2γt)`, which is the convention that
    # lets it be compared with `γ_cold` directly. Drawn with one γ the dashed
    # line peeled a decade below the curve it is supposed to lie along by the
    # end of the window -- 6.1x, 11.6x and 18.0x at these three wavenumbers --
    # so the plot showed a correct fit failing.
    i0 = argmin(abs.(t .- t0))
    tt = range(t0, t1; length = 50)
    plot!(growth, tt, ε_e[i0] .* exp.(2 .* γ .* (tt .- t0));
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
