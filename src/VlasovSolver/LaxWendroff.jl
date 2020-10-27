module LaxWendroff

function generate_solver(f₀, f, c)
    function lax_wendroff!(g, f, i, i⁻, i⁺)
        g[i] = f[i] - 0.5*c*(f[i⁺]-f[i⁻]) + 0.5*c^2*(f[i⁺]-2*f[i]+f[i⁻])
    end

    function solve!()
        lax_wendroff!(f, f₀, 1, length(f), 2)
        @inbounds @fastmath @simd for i = 2:length(f)-1
            lax_wendroff!(f, f₀, i, i-1, i+1)
        end
        lax_wendroff!(f, f₀, length(f), length(f)-1, 1)
    end

    function solve!(h, h₀)
        lax_wendroff!(h, h₀, 1, length(f), 2)
        @inbounds @fastmath @simd for i = 2:length(f)-1
            lax_wendroff!(h, h₀, i, i-1, i+1)
        end
        lax_wendroff!(h, h₀, length(f), length(f)-1, 1)
    end

    return solve!
end

function generate_solver(f₀, f)
    function lax_wendroff!(g, f, i, i⁻, i⁺, c)
        g[i] = f[i] - 0.5*c*(f[i⁺]-f[i⁻]) + 0.5*c^2*(f[i⁺]-2*f[i]+f[i⁻])
    end

    function solve!(c)
        lax_wendroff!(f, f₀, 1, length(f), 2, c)
        @inbounds @fastmath @simd for i = 2:length(f)-1
            lax_wendroff!(f, f₀, i, i-1, i+1, c)
        end
        lax_wendroff!(f, f₀, length(f), length(f)-1, 1, c)
    end

    function solve!(h, h₀, c)
        lax_wendroff!(h, h₀, 1, length(f), 2, c)
        @inbounds @fastmath @simd for i = 2:length(f)-1
            lax_wendroff!(h, h₀, i, i-1, i+1, c)
        end
        lax_wendroff!(h, h₀, length(f), length(f)-1, 1, c)
    end

    return solve!
end
end # module
