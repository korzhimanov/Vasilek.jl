include(joinpath("..","..","src","MaxwellSolver","PoissonFourier1D.jl"))
import .PoissonFourier1D

begin
    Δx = 0.01
    ρ = [sin(2π*i*Δx) for i = 0:1000]
    e = similar(ρ)
    PoissonFourier1D.solve_poisson!(e, ρ, Δx)
    print(sum(@. (e - [-cos(2π*i*Δx)/2π for i = 0:1000])^2))
    @test sum(@. (e - [-cos(2π*i*Δx)/2π for i = 0:1000])^2) ≈ 0 atol=1e-3
end
