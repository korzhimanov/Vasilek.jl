module PFCNonUniform
export make_advect_1D!

function make_advect_1D!(Δx_::Vector{Float64})
    ξ = 2.0
    Δx = copy(Δx_)

    min_Δx = minimum(Δx)/maximum(Δx)
    ξ = (1.0 + min_Δx)*(1.0 + 2*min_Δx)/(3.0 + (min_Δx - 1.0/min_Δx)^2)
    
    f_tmp::Vector{Float64} = similar(Δx)
    
    function ϵ⁺(f::Float64, g::Float64)::Float64
        if f < g
            return min(g - f, ξ*f)
        else
            return max(g - f, -ξ*(1 - f))
        end
    end

    function ϵ⁻(f::Float64, g::Float64)::Float64
        if f < g
            return max(f - g, -ξ*f)
        else
            return min(f - g, ξ*(1 - f))
        end
    end
            
    Φ⁺(α, f, Δx) = α*(f[2] + (Δx[2] - α)/(Δx[3] + Δx[2] + Δx[1])*(
                    ϵ⁺(f[2], f[3])/(Δx[3] + Δx[2])*(Δx[2] + Δx[1] - α) +
                    ϵ⁻(f[2], f[1])/(Δx[2] + Δx[1])*(Δx[3] + α)))
    
    Φ⁻(α, f, Δx) = α*(f[2] - (Δx[2] + α)/(Δx[3] + Δx[2] + Δx[1])*(
                    ϵ⁺(f[2], f[3])/(Δx[3] + Δx[2])*(Δx[1] - α) +
                    ϵ⁻(f[2], f[1])/(Δx[2] + Δx[1])*(Δx[2] + Δx[3] + α)))
    
    function advect_1D!(f, α::Float64)
        f_tmp[:] = copy(f)
        if α > 0
            Φ = Φ⁺(α, [f[end],f[1],f[2]], [Δx[end],Δx[1],Δx[2]])
            f_tmp[1] = f_tmp[1] - Φ/Δx[1]
            f_tmp[2] = f_tmp[2] + Φ/Δx[2]
            for i in 2:length(Δx)-1
                Φ = Φ⁺(α, f[i-1:i+1], Δx[i-1:i+1])
                f_tmp[i] = f_tmp[i] - Φ/Δx[i]
                f_tmp[i+1] = f_tmp[i+1] + Φ/Δx[i+1]
            end
            Φ = Φ⁺(α, [f[end-1],f[end],f[1]], [Δx[end-1],Δx[end],Δx[1]])
            f_tmp[end] = f_tmp[end] - Φ/Δx[end]
            f_tmp[1] = f_tmp[1] + Φ/Δx[1]
        else
            Φ = Φ⁻(α, [f[end],f[1],f[2]], [Δx[end],Δx[1],Δx[2]])
            f_tmp[1] = f_tmp[1] + Φ/Δx[1]
            f_tmp[end] = f_tmp[end] - Φ/Δx[end]
            for i in 2:length(Δx)-1
                Φ = Φ⁻(α, f[i-1:i+1], Δx[i-1:i+1])
                f_tmp[i] = f_tmp[i] + Φ/Δx[i]
                f_tmp[i-1] = f_tmp[i-1] - Φ/Δx[i-1]
            end
            Φ = Φ⁻(α, [f[end-1],f[end],f[1]], [Δx[end-1],Δx[end],Δx[1]])
            f_tmp[end] = f_tmp[end] + Φ/Δx[end]
            f_tmp[end-1] = f_tmp[end-1] - Φ/Δx[end-1]
        end
        f[:] = copy(f_tmp)
    end
    
    return advect_1D!
end

end
