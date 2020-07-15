module BGK

using NumericalIntegration

function generate_solver(f, v_, Δt, τ)
    e = exp(-Δt/τ)
    v = copy(v_)

    function solve!()
        n = integrate(v, f)
        u = integrate(v, v.*f)./n
        T = integrate(v, @. (v-u)^2 * f)./n

        M = @. n/sqrt(2π*T)*exp(-(v-u)^2/(2T))

        f[:] = f*e + (1.0 - e)*M
    end

    return solve!
end
end  # module BGK
