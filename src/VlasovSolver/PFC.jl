module PFC

function generate_solver(f₀, f)
    fₘₐₓ = maximum(f₀)

    function ϵ⁺(f, g)
        if f < g
            return min(g - f, 2.0*f)
        else
            return max(g - f, -2.0*(fₘₐₓ - f))
        end
    end

    function ϵ⁻(f, g)
        if f < g
            return max(f - g, -2.0*f)
        else
            return min(f - g, 2.0*(fₘₐₓ - f))
        end
    end

    Φ⁺(f, i, i⁻, i⁺, c) = c*(f[i] + (1.0 - c)/3.0*(
                    ϵ⁺(f[i], f[i⁺])/2.0*(2.0 - c) +
                    ϵ⁻(f[i], f[i⁻])/2.0*(1.0 + c)))

    Φ⁻(f, i, i⁻, i⁺, c) = c*(f[i] - (1.0 + c)/3.0*(
                    ϵ⁺(f[i], f[i⁺])/2.0*(1.0 - c) +
                    ϵ⁻(f[i], f[i⁻])/2.0*(2.0 + c)))

    function solve!(c)
        if c > 0
            Φ = Φ⁺(f₀, 1, length(f), 2, c)
            f[1] = f₀[1] - Φ
            f[2] = f₀[2] + Φ
            for i in 2:length(f)-1
                Φ = Φ⁺(f₀, i, i-1, i+1, c)
                f[i] = f[i] - Φ
                f[i+1] = f₀[i+1] + Φ
            end
            Φ = Φ⁺(f₀, length(f), length(f)-1, 1, c)
            f[end] = f[end] - Φ
            f[1] = f[1] + Φ
        else
            Φ = Φ⁻(f₀, 1, length(f), 2, c)
            f[1] = f₀[1] + Φ
            f[end] = f₀[end] - Φ
            for i in 2:length(f)-1
                Φ = Φ⁻(f₀, i, i-1, i+1, c)
                f[i] = f₀[i] + Φ
                f[i-1] = f[i-1] - Φ
            end
            Φ = Φ⁻(f₀, length(f), length(f)-1, 1, c)
            f[end] = f[end] + Φ
            f[end-1] = f[end-1] - Φ
        end
    end

    return solve!
end

function generate_solver(f₀, f, c)
    fₘₐₓ = maximum(f₀)

    function ϵ⁺(f, g)
        if f < g
            return min(g - f, 2.0*f)
        else
            return max(g - f, -2.0*(fₘₐₓ - f))
        end
    end

    function ϵ⁻(f, g)
        if f < g
            return max(f - g, -2.0*f)
        else
            return min(f - g, 2.0*(fₘₐₓ - f))
        end
    end

    Φ⁺(f, i, i⁻, i⁺, c) = c*(f[i] + (1.0 - c)/3.0*(
                    ϵ⁺(f[i], f[i⁺])/2.0*(2.0 - c) +
                    ϵ⁻(f[i], f[i⁻])/2.0*(1.0 + c)))

    Φ⁻(f, i, i⁻, i⁺, c) = c*(f[i] - (1.0 + c)/3.0*(
                    ϵ⁺(f[i], f[i⁺])/2.0*(1.0 - c) +
                    ϵ⁻(f[i], f[i⁻])/2.0*(2.0 + c)))

    function solve⁺!()
        Φ = Φ⁺(f₀, 1, length(f), 2, c)
        f[1] = f₀[1] - Φ
        f[2] = f₀[2] + Φ
        for i in 2:length(f)-1
            Φ = Φ⁺(f₀, i, i-1, i+1, c)
            f[i] = f[i] - Φ
            f[i+1] = f₀[i+1] + Φ
        end
        Φ = Φ⁺(f₀, length(f), length(f)-1, 1, c)
        f[end] = f[end] - Φ
        f[1] = f[1] + Φ
    end

    function solve⁻!()
        Φ = Φ⁻(f₀, 1, length(f), 2, c)
        f[1] = f₀[1] + Φ
        f[end] = f₀[end] - Φ
        for i in 2:length(f)-1
            Φ = Φ⁻(f₀, i, i-1, i+1, c)
            f[i] = f₀[i] + Φ
            f[i-1] = f[i-1] - Φ
        end
        Φ = Φ⁻(f₀, length(f), length(f)-1, 1, c)
        f[end] = f[end] + Φ
        f[end-1] = f[end-1] - Φ
    end

    if c > 0
        return solve⁺!
    else
        return solve⁻!
    end
end

end  # module PFC
