module Godunov

function Φ⁺_constant(f, i, i⁻, c)
    return c*f[i⁻]
end

function Φ⁻_constant(f, i, i⁻, c)
    return c*f[i]
end

function Φ_linear(f, i, i⁻, c)
    return 0.5*c*(f[i]+f[i⁻])
end

function generate_Φ⁺(Riemann_solver_type)
    if Riemann_solver_type == :Riemann_constant
        return Φ⁺_constant
    elseif Riemann_solver_type == :Riemann_linear
        return Φ_linear
    end
end

function generate_Φ⁻(Riemann_solver_type)
    if Riemann_solver_type == :Riemann_constant
        return Φ⁻_constant
    elseif Riemann_solver_type == :Riemann_linear
        return Φ_linear
    end
end

function generate_solver(f₀, f, Riemann_solver_type)

    Φ⁺ = generate_Φ⁺(Riemann_solver_type)
    Φ⁻ = generate_Φ⁻(Riemann_solver_type)

    function godunov⁺!(g, f, i, i⁻, i⁺, c)
        g[i] = f[i] + Φ⁺(f, i, i⁻, c) - Φ⁺(f, i⁺, i, c)
    end

    function godunov⁻!(g, f, i, i⁻, i⁺, c)
        g[i] = f[i] + Φ⁻(f, i, i⁻, c) - Φ⁻(f, i⁺, i, c)
    end

    function solve⁺!(c)
        godunov⁺!(f, f₀, 1, length(f), 2, c)
        @inbounds @fastmath @simd for i = 2:length(f)-1
            godunov⁺!(f, f₀, i, i-1, i+1, c)
        end
        godunov⁺!(f, f₀, length(f), length(f)-1, 1, c)
    end

    function solve⁻!(c)
        godunov⁻!(f, f₀, 1, length(f), 2, c)
        @inbounds @fastmath @simd for i = 2:length(f)-1
            godunov⁻!(f, f₀, i, i-1, i+1, c)
        end
        godunov⁻!(f, f₀, length(f), length(f)-1, 1, c)
    end

    function solve!(c)
        if c > 0
            solve⁺!(c)
        else
            solve⁻!(c)
        end
    end

    return solve!
end

end # module
