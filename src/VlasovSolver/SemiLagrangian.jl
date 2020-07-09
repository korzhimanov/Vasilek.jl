module SemiLagrangian

using Interpolations

function generate_solver(f₀, f; interpolation_order = :Cubic)
    if interpolation_order == :Linear
        interpolation_function = BSpline(Linear())
    elseif interpolation_order == :Quadratic
        interpolation_function = BSpline(Quadratic(Periodic(OnCell())))
    elseif interpolation_order == :Cubic
        interpolation_function = BSpline(Cubic(Periodic(OnCell())))
    end

    function solve!(c)
        if interpolation_order == :Linear
            itp = interpolate([f₀; f₀[1]], interpolation_function)
        else
            itp = interpolate(f₀, interpolation_function)
        end
        etp = extrapolate(itp, Periodic(OnCell()))

        for i = 1:length(f₀)
            if 1.0 ≤ i - c ≤ length(itp)
                f[i] = itp(i - c)
            else
                f[i] = etp(i - c)
            end
        end
    end

    return solve!
end

end  # module SemiLagrangian
