module BGK

using NumericalIntegration

"""
    generate_solver(f, v_, Δt, τ)

BGK relaxation towards the local Maxwellian, in place on `f`.

The moments and the Maxwellian used to be built from fresh temporaries on
every call (`v.*f`, `(v-u)^2 .* f`, `M`, and three more inside
`f*e + (1-e)*M`) — 39 kB per step at 801 velocity points. The two buffers are
now allocated once and every expression is a fused in-place broadcast, which
leaves the arithmetic per element unchanged.
"""
function generate_solver(f, v_, Δt, τ)
    e = exp(-Δt/τ)
    v = copy(v_)

    M = similar(f)
    moment = similar(f)

    function solve!()
        n = integrate(v, f)

        @. moment = v*f
        u = integrate(v, moment)/n

        @. moment = (v - u)^2 * f
        T = integrate(v, moment)/n

        @. M = n/sqrt(2π*T)*exp(-(v - u)^2/(2T))
        @. f = f*e + (1.0 - e)*M
    end

    return solve!
end
end  # module BGK
