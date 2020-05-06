module Upwind

function generate_solver(f₀, f, c)
   function upwind⁺!(g, f, i, i⁻)
        g[i] = f[i] - c*(f[i]-f[i⁻])
    end

    function solve⁺!()
        upwind⁺!(f, f₀, 1, length(f))
        @inbounds @fastmath @simd for i = 2:length(f)
            upwind⁺!(f, f₀, i, i-1)
        end
    end

    function upwind⁻!(g, f, i, i⁺)
        g[i] = f[i] - c*(f[i⁺]-f[i])
    end

    function solve⁻!()
        @inbounds @fastmath @simd for i = 1:length(f)-1
            upwind⁻!(f, f₀, i, i+1)
        end
        upwind⁻!(f, f₀, length(f), 1)
    end

    if c > 0
        return solve⁺!
    else
        return solve⁻!
    end
end

function generate_solver(f₀, f)
   function upwind⁺!(g, f, i, i⁻, c)
        g[i] = f[i] - c*(f[i]-f[i⁻])
    end

    function upwind⁻!(g, f, i, i⁺, c)
        g[i] = f[i] - c*(f[i⁺]-f[i])
    end

    function solve⁺!(c)
        upwind⁺!(f, f₀, 1, length(f), c)
        @inbounds @fastmath @simd for i = 2:length(f)
            upwind⁺!(f, f₀, i, i-1, c)
        end
    end

    function solve⁻!(c)
        @inbounds @fastmath @simd for i = 1:length(f)-1
            upwind⁻!(f, f₀, i, i+1, c)
        end
        upwind⁻!(f, f₀, length(f), 1, c)
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
