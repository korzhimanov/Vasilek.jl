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
    
    function advect_1D!(f, α::Float64)
        f_tmp[:] = copy(f)
        if α > 0
            for i in 2:length(Δx)-1
                Φ = α*(f[i] + (Δx[i] - α)/(Δx[i+1] + Δx[i] + Δx[i-1])*(
                    ϵ⁺(f[i], f[i+1])/(Δx[i+1] + Δx[i])*(Δx[i] + Δx[i-1] - α) +
                    ϵ⁻(f[i], f[i-1])/(Δx[i] + Δx[i-1])*(Δx[i+1] + α)))
                f_tmp[i] = f_tmp[i] - Φ/Δx[i]
                f_tmp[i+1] = f_tmp[i+1] + Φ/Δx[i+1]
            end
        else
            for i in 2:length(Δx)-1
                Φ = α*(f[i] - (Δx[i] + α)/(Δx[i+1] + Δx[i] + Δx[i-1])*(
                    ϵ⁺(f[i], f[i+1])/(Δx[i+1] + Δx[i])*(Δx[i-1] - α) +
                    ϵ⁻(f[i], f[i-1])/(Δx[i] + Δx[i-1])*(Δx[i] + Δx[i+1] + α)))
                f_tmp[i] = f_tmp[i] + Φ/Δx[i]
                f_tmp[i-1] = f_tmp[i-1] - Φ/Δx[i-1]
            end
        end
        f[:] = copy(f_tmp)
    end
    
    return advect_1D!
end

end
