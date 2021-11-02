module FDTD1D
export make_advance_fields, YeeMesh1D

struct YeeMesh1D{T<:AbstractFloat, S<:Integer}
    ey::Vector{T}
    ez::Vector{T}
    hy::Vector{T}
    hz::Vector{T}
    N::S
    function YeeMesh1D{T}(N::S) where {T, S}
        new{T,S}(zeros(T, N+1), zeros(T, N+1), zeros(T, N), zeros(T, N), N)
    end
end

function make_advance_fields(f::YeeMesh1D{T,S}, cfl, pulse_shape, Δt, Δx, x_min) where {T,S}
    Nx = f.N

    function generate_fields_x_min!(t)
        f.ey[2] -= pulse_shape["y"].(t, x_min + Δx)
        f.ez[2] -= pulse_shape["z"].(t, x_min + Δx)
        
        f.hz[2] -= pulse_shape["y"].(t + 0.5*Δt, x_min + 1.5*Δx)
        f.hy[2] -= pulse_shape["z"].(t + 0.5*Δt, x_min + 1.5*Δx)
    end

    function update_ey!()
        for i = 2:Nx
            f.ey[i] -= cfl*(f.hz[i] - f.hz[i-1])
        end
    end
    
    function update_ez!()
        for i = 2:Nx
            f.ez[i] += cfl*(f.hy[i] - f.hy[i-1])
        end
    end
    
    function update_hy!()
        for i = 1:Nx
            f.hy[i] += cfl*(f.ez[i+1] - f.ez[i])
        end
    end
    
    function update_hz!()
        for i = 1:Nx
            f.hz[i] -= cfl*(f.ey[i+1] - f.ey[i])
        end
    end

    function make_step!(j)
        update_ey!()
        update_ez!()
        
        update_hz!()
        update_hy!()
    end
    
    function advance_fields!(t, j)
        generate_fields_x_min!(t)
        make_step!(j)
    end
    
    return advance_fields!
end

end
