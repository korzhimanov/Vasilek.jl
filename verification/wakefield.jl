using Plots
using NumericalIntegration

include(joinpath("..", "src","VlasovSolver","PFC.jl"))
include(joinpath("..", "src","MaxwellSolver","FDTD1D.jl"))
include(joinpath("..", "src","MaxwellSolver","PoissonFourier1D.jl"))

Δx = 0.1*2π
Δt = 0.05*Δx
Δp = 0.1

x_min = -5.0*2π
box_length = 20.0*2π
plasma_thickness = 10.0*2π
plasma_temperature = 0.1
plasma_density = 0.1

laser_amplitude = 1.0
laser_duration = 5*2π

total_time = 2π*10

x = collect(x_min:Δx:(box_length+x_min))
p = collect(-laser_amplitude*4:Δp:laser_amplitude*4)

Nx = length(x)
Np = length(p)

f = plasma_density/sqrt(2π*plasma_temperature) *
    (@. exp(-0.5*(p)^2/plasma_temperature)) *
    (@. 0.5*(tanh(x)-tanh(x-plasma_thickness)))'
    #  (@. (x+0.1)./(x+0.1))'
ni = integrate(p, f)
Ni = integrate(x, ni)

n0 = integrate(p, f)
Ne = integrate(x, n0)

g = copy(f)

advect_x! = PFC.generate_solver(g[1,:], f[1,:]; fₘₐₓ = maximum(f))
advect_p! = PFC.generate_solver(f[:,1], g[:,1]; fₘₐₓ = maximum(f))

em = FDTD1D.YeeMesh1D{Float64}(length(x)-1)

pulse_shape = (
    y = (t,x) -> exp(-((x-t)/laser_duration)^2)*sin(x-t),
    z = (t,x) -> 0.0
)

advance_fields! = FDTD1D.make_advance_fields(em, Δt/Δx, pulse_shape, Δt, Δx, x_min, FDTD1D.PML(0,1.0,Δx,Δt))

t = collect(0.0:Δt:total_time)
n = t * n0'
n[1,:] = n0

ε_e = similar(t)
ε = similar(t)

solve_poisson! = PoissonFourier1D.generate_solver(n0, Δx)

e = similar(x)
solve_poisson!(e, n[1,:]-ni)

ε_e[1] = integrate(x, e.^2)
ε[1] = integrate(x, integrate(p, @. f*p^2)) + ε_e[1]

py = zeros(length(x))
pz = zeros(length(x))

for k in 2:length(t)
    
    for j = 1:Np
        advect_x!(view(g,j,:), view(f,j,:), p[j]*Δt/Δx)
    end

    n[k,:] = integrate(p, g)
    solve_poisson!(e, n[k,:]-ni)

    py[:] = @. py + em.ey*Δt
    pz[:] = @. pz + em.ez*Δt

    advance_fields!(k*Δt, (y=py.*n[k,:], z=pz.*n[k,:]))

    for i = 1:Nx
        advect_p!(view(f,:,i), view(g,:,i), e[i]*Δt/Δp)
    end

    ε_e[k] = integrate(x, e.^2)
    ε[k] = integrate(x, integrate(p, @. f*p^2)) + ε_e[k]
end

heatmap(x/2π, t/2π, em.ey)