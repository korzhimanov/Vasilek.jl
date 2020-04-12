module LaxWendroff

function generate_solver(f₀, f, c)
    function lax_wendroff!(g, f, i, i⁻, i⁺)
        g[i] = f[i] - 0.5*c*(f[i⁺]-f[i⁻]) + 0.5*c^2*(f[i⁺]-2*f[i]+f[i⁻])
    end

    function solve!()
        lax_wendroff!(f, f₀, 1, length(f), 2)
        for i = 2:length(f)-1
            # g[i] = f[i] - c*(f[i+1]-f[i-1]) + 2*c^2*(f[i+1]-2*f[i]+f[i-1])
            lax_wendroff!(f, f₀, i, i-1, i+1)
        end
        lax_wendroff!(f, f₀, length(f), length(f)-1, 1)
    end
    return solve!
end

end
