# The Vlasov–Poisson driver ships with the package. These run without the
# verification harness, as a user would: nothing here is included from test/.

using Vasilek.VlasovPoisson1D1V: line_advector, cell_widths, substeps

@testset "vlasov_poisson is the package's" begin
    @test isdefined(Vasilek, :vlasov_poisson)
    @test :vlasov_poisson in names(Vasilek)
    @test Vasilek.vlasov_poisson === Vasilek.VlasovPoisson1D1V.vlasov_poisson

    # A short Landau run: the field energy decays, the mass and the total
    # energy hold, and the requested mode is recorded.
    k = 0.5
    Nx = 32
    x = collect(range(2π/k/Nx; step = 2π/k/Nx, length = Nx))
    v = collect(-4.0:0.2:4.0)
    f₀ = [exp(-u^2/2)/sqrt(2π)*(1 + 0.01*cos(k*y)) for u in v, y in x]
    t = collect(0.0:0.1:10.0)
    r = vlasov_poisson(x, v, f₀, t; modes = (k,), invariants = true)
    @test size(r.f) == size(f₀)
    @test size(r.E_modes) == (length(t), 1)
    @test r.ε_e[end] < r.ε_e[1]
    @test maximum(abs, r.mass .- r.mass[1])/r.mass[1] < 1e-12
    @test maximum(abs, r.ε .- r.ε[1])/r.ε[1] < 1e-3
    @test minimum(r.fmin) ≥ 0

    # the adapter that takes a displacement whatever the scheme takes
    col = [1.0 + 0.5sin(2π*i/32) for i in 1:32]
    ref = advect!(similar(col), col, Upwind(), 0.25)
    @test line_advector(Upwind(), fill(2.0, 32))(copy(col), 0.5) == ref
    @test_throws ErrorException line_advector(Upwind(), [1.0, 2.0, 1.0, 1.0])
    @test substeps(2.5, 1.0) == 3
    @test cell_widths([0.0, 1.0, 3.0]) == [1.0, 1.5, 2.0]

    # The field solve needs a uniform x grid. A stretched one used to run and
    # return a field 16% off in energy without a word; it is now refused, and
    # a collected range, whose spacing rounds, is still taken.
    s = (0:Nx-1)/Nx
    stretched = @. 2π/k*(s + 0.08sin(2π*s)/(2π)) + x[1]
    @test_throws ArgumentError vlasov_poisson(stretched, v, f₀, t)
    @test_throws ArgumentError Vasilek.VlasovPoisson1D1V.make_poisson(stretched)
    @test Vasilek.VlasovPoisson1D1V.make_poisson(collect(range(0.1, 25.0, length = 97))) isa Function
    # Uniform is to rounding, which is on the scale of the largest coordinate
    # rather than of the spacing. A Float32 grid resolves its spacing to a few
    # parts in 10⁶, which 1e-10 of the spacing alone refuses although the
    # harness's driver ran it; it runs, and matches the Float64 run to Float32
    # precision (measured 5.6e-7 of the peak in `f`, 1.8e-6 in `ε_e`). A grid
    # far from the origin is taken too, and a stretch of a millionth of a cell
    # is still refused.
    r₃₂ = vlasov_poisson(Float32.(x), Float32.(v), Float32.(f₀), Float32.(t))
    @test maximum(abs, r₃₂.f .- r.f) < 1e-5*maximum(r.f)
    @test maximum(abs, r₃₂.ε_e .- r.ε_e) < 1e-4*maximum(r.ε_e)
    @test Vasilek.VlasovPoisson1D1V.make_poisson(collect(1e6 .+ (0:63)*0.1)) isa Function
    @test_throws ArgumentError Vasilek.VlasovPoisson1D1V.make_poisson(@. x + 1e-6*(2π/k/Nx)*sin(2π*s))
    # A single time has no step to record a history in: refused, rather than a
    # BoundsError from copying the last entry forward.
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, [0.0])
end
