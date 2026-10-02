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

# A scheme the driver has never heard of. As a partner for a default it is
# refused whatever it would do, since nothing vouches for its bounds; it never
# takes a step here, so it needs no `advect!`.
struct UnvouchedScheme <: AbstractAdvection1D end

using Vasilek.VlasovPoisson1D1V: keeps_bounds

@testset "a scheme given alone has to keep the default's bounds" begin
    # One wavelength of k = 0.5 at 50% amplitude, 64 × 121 over ±6. With the
    # other direction left to its default, a PFCNonUniform on [0, maximum(f)],
    # a non-positive scheme used to stop the run from inside the default with a
    # DomainError: LaxWendroff in v handed it f = -1.9e-9 on the 16th step, the
    # cubic SemiLagrangian -7.0e-10 on the 21st. Both calls are now an
    # ArgumentError before the first step -- in a one-step run too, which the old
    # failure could not reach -- in either direction, and with collisions.
    x = collect(range(4π/64; step = 4π/64, length = 64))
    v = collect(-6.0:0.1:6.0)
    f₀ = [exp(-u^2/2)/sqrt(2π)*(1 + 0.5cos(0.5y)) for u in v, y in x]
    t = collect(0.0:0.05:50.0)
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, t; scheme_v = LaxWendroff())
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, t; scheme_v = SemiLagrangian(CubicSpline()))
    for s in (LaxWendroff(), SemiLagrangian(CubicSpline()))
        @test_throws ArgumentError vlasov_poisson(x, v, f₀, t[1:2]; scheme_v = s)
        @test_throws ArgumentError vlasov_poisson(x, v, f₀, t[1:2]; scheme_x = s)
        @test_throws ArgumentError vlasov_poisson(x, v, f₀, t[1:2]; scheme_v = s,
                                                  collisions = BGK(1.0))
    end
    message = try
        vlasov_poisson(x, v, f₀, t[1:2]; scheme_v = LaxWendroff())
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("LaxWendroff()", message) && occursin("Pass scheme_x as well", message)
    # Refused as well: a PFC bounded wider than [0, maximum(f)], whose limiter
    # lets f out to its own bound, and a scheme of a type the driver does not know.
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, t[1:2]; scheme_v = PFC(fmin = 0.0, fmax = 1.0))
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, t[1:2]; scheme_x = UnvouchedScheme())

    # Given both schemes, the driver checks neither. The default given
    # explicitly runs as it always did, and stops as the refusal says it would;
    # the partners the message suggests run on, keep the mass, and leave f below
    # 0 where LaxWendroff takes it.
    Δx = cell_widths(x)
    bounded = f -> PFCNonUniform(Δx; fmin = 0.0, fmax = maximum(f))
    @test_throws DomainError vlasov_poisson(x, v, f₀, t[1:41]; scheme_v = LaxWendroff(),
                                            scheme_x = bounded)
    for sx in (LaxWendroff(), PFCNonUniform(Δx; fmin = -Inf, fmax = Inf))
        r = vlasov_poisson(x, v, f₀, t[1:41]; scheme_v = LaxWendroff(), scheme_x = sx,
                           invariants = true)
        @test minimum(r.fmin) < 0
        @test maximum(abs, r.mass .- r.mass[1])/r.mass[1] < 1e-12
    end
    # A scheme that keeps the bounds still partners a default, as before.
    for s in (Godunov(PiecewiseLinear(), VanLeer()), f -> PFC(fmin = 0.0, fmax = maximum(f)))
        r = vlasov_poisson(x, v, f₀, t[1:41]; scheme_v = s, invariants = true)
        @test minimum(r.fmin) ≥ 0
    end

    # The classification, and the property it rests on. From a square pulse,
    # twenty steps either way, the schemes it accepts stay inside [0, 1] --
    # exactly, on this pulse; the bound leaves room for round-off -- and the
    # ones it refuses leave it by 0.16 (the cubic spline) to 0.24 (LaxWendroff).
    accepted = (Upwind(), Godunov(PiecewiseConstant()), Godunov(PiecewiseLinear(), VanLeer()),
                Godunov(PiecewiseLinear(), Superbee()), SemiLagrangian(LinearSpline()),
                PFC(fmin = 0.0, fmax = 1.0))
    refused = (LaxWendroff(), Godunov(PiecewiseLinear()), SemiLagrangian(QuadraticSpline()),
               SemiLagrangian(CubicSpline()))
    pulse = [16 < i ≤ 40 ? 1.0 : 0.0 for i in 1:64]
    function excursion(s)
        worst = 0.0
        for c in (-0.77, 0.3)
            src, dst, ws = copy(pulse), similar(pulse), workspace(s, length(pulse))
            for _ in 1:20
                advect!(dst, src, s, c, ws)
                copyto!(src, dst)
                worst = max(worst, -minimum(src), maximum(src) - 1)
            end
        end
        return worst
    end
    for s in accepted
        @test keeps_bounds(s, 0.0, 1.0)
        @test excursion(s) ≤ 4eps()
    end
    for s in refused
        @test !keeps_bounds(s, 0.0, 1.0)
        @test excursion(s) > 0.05
    end
    @test keeps_bounds(PFCNonUniform(fill(0.5, 8); fmin = 0.0, fmax = 1.0), 0.0, 1.0)
    @test !keeps_bounds(PFC(fmin = 0.0, fmax = 1.5), 0.0, 1.0)
    @test !keeps_bounds(PFC(fmin = -0.1, fmax = 1.0), 0.0, 1.0)
    @test !keeps_bounds(UnvouchedScheme(), 0.0, 1.0)
end

@testset "with collisions the defaults' upper bound is lifted" begin
    # BGK relaxes a line towards a Maxwellian whose peak can sit above the
    # line's, so the maximum of f is not kept, and the defaults' upper bound,
    # that maximum, stopped valid runs with a DomainError from inside the step.
    # With `collisions` the defaults are bounded below only.
    x = collect(range(4π/64; step = 4π/64, length = 64))
    v = collect(-6.0:0.1:6.0)
    # A line flat over |v| < 2: under BGK(0.1) it reached 0.416 against 0.375.
    flat = [abs(u) < 2 ? 0.25*(1 + 0.5cos(0.5y)) : 0.0 for u in v, y in x]
    r = vlasov_poisson(x, v, flat, collect(0.0:0.05:5.0); collisions = BGK(0.1),
                       invariants = true)
    @test maximum(r.fmax[1:end-1]) > 1.1*maximum(flat)
    @test minimum(r.fmin) ≥ 0
    @test maximum(abs, r.mass .- r.mass[1])/r.mass[1] < 1e-12
    # A Superbee v sweep, which partners a default, under BGK(1.0) at 50%: it
    # went 3.3e-7 of the maximum above it on the way.
    f₀ = [exp(-u^2/2)/sqrt(2π)*(1 + 0.5cos(0.5y)) for u in v, y in x]
    r = vlasov_poisson(x, v, f₀, collect(0.0:0.05:50.0);
                       scheme_v = Godunov(PiecewiseLinear(), Superbee()),
                       collisions = BGK(1.0), invariants = true)
    @test minimum(r.fmin) ≥ 0
    # The lower bound still holds the partner to it: LaxWendroff given alone is
    # refused with collisions as without, and a PFC bounded wider above is now
    # taken, since nothing above is bounded.
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, [0.0, 0.05]; scheme_v = LaxWendroff(),
                                              collisions = BGK(1.0))
    @test vlasov_poisson(x, v, f₀, [0.0, 0.05]; scheme_v = PFC(fmin = 0.0, fmax = 1.0),
                         collisions = BGK(1.0)).f isa Matrix
    # Without collisions nothing changes: the same PFC is refused.
    @test_throws ArgumentError vlasov_poisson(x, v, f₀, [0.0, 0.05]; scheme_v = PFC(fmin = 0.0, fmax = 1.0))
end
