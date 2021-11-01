module FDTD1D
export make_advance_fields

function make_advance_fields(Nx, half_cfl, pulse_shape, Δt, Δx, x_min)
    
    function generate_fields_x_min!(f, t)
        f["ey"][2] -= pulse_shape["y"].(t, x_min + Δx)
        f["ez"][2] -= pulse_shape["z"].(t, x_min + Δx)
        
        f["hz"][2] -= pulse_shape["y"].(t + 0.5*Δt, x_min + 1.5*Δx)
        f["hy"][2] -= pulse_shape["z"].(t + 0.5*Δt, x_min + 1.5*Δx)
    end

    function update_ey!(ey, hz)
        for i = 2:Nx
            ey[i] -= half_cfl*(hz[i] - hz[i-1])
        end
    end
    
    function update_ez!(ez, hy)
        for i = 2:Nx
            ez[i] += half_cfl*(hy[i] - hy[i-1])
        end
    end
    
    function update_hy!(hy, ez)
        for i = 1:Nx
            hy[i] += half_cfl*(ez[i+1] - ez[i])
        end
    end
    
    function update_hz!(hz, ey)
        for i = 1:Nx
            hz[i] -= half_cfl*(ey[i+1] - ey[i])
        end
    end

    function make_step!(f, j)
        update_ey!(f["ey"], f["hz"])
        update_ez!(f["ez"], f["hy"])
        
        update_hz!(f["hz"], f["ey"])
        update_hy!(f["hy"], f["ez"])
    end
    
    function advance_fields!(f, t, j)
        generate_fields_x_min!(f, t)
        make_step!(f, j)
    end
    
    return advance_fields!
end

end
