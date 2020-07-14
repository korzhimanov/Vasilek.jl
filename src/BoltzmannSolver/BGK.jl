module BGK

using NumericalIntegration

function generate_solver(f, v_, νΔt)
    v = copy(v_)

    function solve!()
        n = integrate(v, f)
        u = integrate(v, v.*f)./n
        T = integrate(v, @. (v-u)^2 * f)./n

        M = @. n/sqrt(2π*T)*exp(-(v-u)^2/(2T))

        f[:] = f*(1 - νΔt) + νΔt*M
    end

    return solve!
end
end  # module BGK
