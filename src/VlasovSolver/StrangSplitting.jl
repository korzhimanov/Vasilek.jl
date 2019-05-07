module StrangSplitting
export make_time_step_2d!

function make_time_step_2d!(f, vΔt, advect!)
    for (fi, α) in zip(eachcol(f[1]), vΔt[1](f[2]))
        advect![1](fi, α/2)
    end
    f[2][:] = (f[1])'
    for (fi, α) in zip(eachcol(f[2]), vΔt[2](f[1]))
        advect![2](fi, α)
    end
    f[1][:] = (f[2])'
    for (fi, α) in zip(eachcol(f[1]), vΔt[1](f[2]))
        advect![1](fi, α/2)
    end
end

end
