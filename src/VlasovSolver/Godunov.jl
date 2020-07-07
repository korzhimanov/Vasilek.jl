module Godunov

include(joinpath(dirname(@__FILE__),"Limiters.jl"))

function calc_limiter(limiter_func, f, i, i⁻, i⁺)
    if f[i] - f[i⁻] ≈ 0.0
        r = 0.0
    elseif f[i⁺] - f[i] ≈ 0.0
        r = 2.0
    else
        r = (f[i] - f[i⁻])/(f[i⁺] - f[i])
    end
    return limiter_func(r)
end

function Φ_constant(f, i, i⁻, i⁻², c)
    return abs(c)*f[i⁻]
end

function generate_Φ_linear(limiter_func)
    function Φ(f, i, i⁻, i⁻², c)
        return abs(c)*(f[i⁻] + calc_limiter(limiter_func, f, i⁻, i⁻², i)*0.5*(f[i]-f[i⁻]))
    end
    return Φ
end

function generate_Φ(Riemann_solver_type, limiter_func)
    if Riemann_solver_type == :Riemann_constant
        return Φ_constant
    elseif Riemann_solver_type == :Riemann_linear
        return generate_Φ_linear(limiter_func)
    end
end

function generate_solve_pm(f₀, f, Riemann_solver_type; flux_limiter)
    if flux_limiter == :None
        limiter_func = x -> 1.0
    else
        limiter_func = Limiters.generate_limiter(flux_limiter)
    end
    Φ = generate_Φ(Riemann_solver_type, limiter_func)

    function godunov⁺!(g, f, i, i⁻, i⁺, i⁻², c)
        g[i] = f[i] + Φ(f, i, i⁻, i⁻², c) - Φ(f, i⁺, i, i⁻, c)
    end

    function godunov⁻!(g, f, i, i⁻, i⁺, i⁺², c)
        return godunov⁺!(g, f, i, i⁺, i⁻, i⁺², c)
    end

    function solve⁺!(c)
        godunov⁺!(f, f₀, 1, length(f), 2, length(f)-1, c)
        godunov⁺!(f, f₀, 2, 1, 3, length(f), c)
        @inbounds @fastmath @simd for i = 3:length(f)-1
            godunov⁺!(f, f₀, i, i-1, i+1, i-2, c)
        end
        godunov⁺!(f, f₀, length(f), length(f)-1, 1, length(f)-2, c)
    end

    function solve⁻!(c)
        godunov⁻!(f, f₀, 1, length(f), 2, 3, c)
        @inbounds @fastmath @simd for i = 2:length(f)-2
            godunov⁻!(f, f₀, i, i-1, i+1, i+2, c)
        end
        godunov⁻!(f, f₀, length(f)-1, length(f)-2, length(f), 1, c)
        godunov⁻!(f, f₀, length(f), length(f)-1, 1, 2, c)
    end

    return solve⁻!, solve⁺!
end

function generate_solver(f₀, f, Riemann_solver_type; flux_limiter=:None)
    solve⁻!, solve⁺! = generate_solve_pm(f₀, f, Riemann_solver_type; flux_limiter=flux_limiter)

    function solve!(c)
        if c > 0
            solve⁺!(c)
        else
            solve⁻!(c)
        end
    end

    return solve!
end

function generate_solver(f₀, f, c, Riemann_solver_type; flux_limiter=:None)
    solve⁻!, solve⁺! = generate_solve_pm(f₀, f, Riemann_solver_type; flux_limiter=flux_limiter)

    if c > 0
        return () -> solve⁺!(c)
    else
        return () -> solve⁻!(c)
    end
end

end # module
