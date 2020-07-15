MODULES  = (:BGK,
            :Landau1P)

for mod_name in MODULES
    include(joinpath(dirname(@__FILE__),"..","..","src","BoltzmannSolver","$mod_name.jl"))
end

function test_1v_damping_step(mod_name, Δt, f₀, f₁, v, exp_norm_dev, args...; plot_needed=false, kwargs...)
    if plot_needed
        plt = plot(v, f₀, label = "initial")
        plot!(v, f₁, label = "expected")
    end

    f = copy(f₀)

    if mod_name == :BGK
        damp! = eval(:($mod_name)).generate_solver(f, v, Δt, args...; kwargs...)
    else
        damp! = eval(:($mod_name)).generate_solver(f₀, f, v, Δt, args...; kwargs...)
    end
    for k = 1:100
        damp!()
        # plot!(v, f)
        f₀ = copy(f)
    end

    if plot_needed
        plot!(v, f, label = "$mod_name $args")
    end

    println("$mod_name $args $(values(values(kwargs))): $(norm(f - f₁))")
    @test norm(f - f₁) ≈ 0 atol=exp_norm_dev

    if plot_needed
        display(plt)
    end
end

@testset "Test 1V Boltzmann solvers" begin
    Δt = 0.1
    v = collect(-4:0.1:4)

    f₀ = @. exp(-v^2)
    f₁ = @. exp(-v^2)

    test_1v_damping_step(:BGK, Δt, f₀, f₁, v, 1e-4, 1e-2)
    test_1v_damping_step(:Landau1P, Δt, f₀, f₁, v, 3e-3, 1e-2)

    a = 3e-1
    n = 0.5*sqrt(π)*(a+2.)
    T = 0.25*sqrt(π)*(3a+2.)/n
    f₀ = @. exp(-v^2)*(1.0 + a*v^2)
    f₁ = @. n/sqrt(2π*T)*exp(-v^2/(2T))

    test_1v_damping_step(:BGK, Δt, f₀, f₁, v, 2e-3, 1e-1)
end
