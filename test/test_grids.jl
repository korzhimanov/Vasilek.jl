using Vasilek
using Vasilek.Grids
using Vasilek.Advection: substep_bound, nsubsteps, takes, Courant, Displacement

# A scheme of one's own that takes a displacement: `takes` is the extension point.
struct DisplacedShift <: AbstractAdvection1D end
Vasilek.Advection.takes(::DisplacedShift) = Displacement()

# Grids are types: whether a grid is uniform is decided where it is built, by
# `grid` and the constructors, and everything after dispatches on the answer.
# These are their rules; the schemes on them are tested beside `OnGrid` in
# `test/VlasovSolver/test_strang_splitting.jl` and the driver in `test_driver.jl`.
@testset "Grids" begin
    @testset "a range is a uniform grid, its step exact" begin
        for r in (-4.0:0.2:4.0, range(-6, 6; length = 41), 1e6 .+ (0:0.1:10),
                  range(4π/32; step = 4π/32, length = 32), -4f0:0.1f0:4f0,
                  -4//1:1//5:4//1, 1:10)
            T = float(eltype(r))
            g = grid(r)
            @test g isa UniformGrid{T} && UniformGrid(r) isa UniformGrid{T}
            @test spacing(g) === T(step(r))
            @test nodes(g) == collect(T, r) && length(g) == length(r)
        end
        # The collected nodes' first difference is not the step the range was
        # written with; the range carries it to the bit, and a vector's spacing
        # is its mean, `(z[end] - z[1])/(n - 1)`, which here is the step too.
        @test spacing(grid(-4.0:0.2:4.0)) === 0.2
        let v = collect(-4.0:0.2:4.0)
            @test v[2] - v[1] === 0.20000000000000018 && spacing(grid(v)) === 0.2
        end
        # A range's cell widths are its step in every cell, as its grid's are,
        # not the differences of its nodes, which round in the last bit.
        @test cell_widths(-4.0:0.2:4.0) == fill(0.2, 41) == widths(grid(-4.0:0.2:4.0))
        @test cell_widths(collect(-4.0:0.2:4.0)) != fill(0.2, 41)
        @test_throws ArgumentError grid(10.0:-1.0:1.0)     # nodes have to increase
        # A range's nodes are where the run samples its data, and are held to
        # the cap as a vector's are: a Float16 range over [−4, 4] rounds its
        # nodes by 1% of a cell, and is not uniform, whatever its step.
        r16 = range(Float16(-4), Float16(4); length = 41)
        @test grid(r16) isa NonUniformGrid{Float16}
        @test_throws ArgumentError UniformGrid(r16)
        @test grid(Float16(-4):Float16(0.25):Float16(4)) isa UniformGrid{Float16}
        # A range's cell widths are its grid's: its step in every cell, or its
        # nodes' own where it is not uniform, which threw.
        @test cell_widths(-4.0:0.2:4.0) == fill(0.2, 41)
        @test cell_widths(r16) == cell_widths(collect(r16))
        @test isbitstype(typeof(UniformGrid(0.1, 8)))
    end

    @testset "a vector of nodes is classified once" begin
        # `grid`'s rule: every spacing the first to 1e-10 of it, or to 8 ulps of
        # the largest node in the nodes' own type, and never past 0.5% of a
        # cell. Pinned here are rows 1.9 times or more inside or outside the
        # bound on Julia 1.13, max|Δ|/h in the comments; the table is in
        # docs/src/driver-notes.md, "Grids".
        L = 4π
        σ = range(-1, 1; length = 41)
        s = (0:31)/32
        # `Any[]`: a literal of vectors promotes them to one type, Float32 nodes
        # to Float64 ones, which are not uniform to Float64 rounding.
        uniform = Any[
            collect(range(L/32; step = L/32, length = 32)),                    # 3.4e-15
            collect(range(L/10^6; step = L/10^6, length = 10^6)),              # 1.4e-10
            collect(LinRange(0.0, L, 10^6)),                                    # 1.8e-10
            cumsum(fill(L/10^6, 10^6)),                                         # 1.4e-10
            collect(range(-Float32(L/2), Float32(L/2); length = 2000)),         # 7.6e-5
            cumsum(fill(Float32(L/2000), 2000)),                                # 2.1e-4
            collect(range(Float32(L/20000); step = Float32(L/20000), length = 20000)),  # 1.3e-3
            Float32.(collect(range(π/8; step = π/8, length = 64))),
            collect(1e6 .+ (0:63)*0.1),                                         # 1.2e-9
            [(j - 0.5)*4π/6 for j in 1:6],                                      # 4.2e-16
            collect(-4//1:1//5:4//1), collect(1:32),                            # exact
            Float16.(collect(-4.0:0.25:4.0)),                                   # exact in Float16
        ]
        nonuniform = Any[
            collect(range(Float32(L/10^6); step = Float32(L/10^6), length = 10^6)),  # 6.3e-2
            Float16.(collect(-4.0:0.2:4.0)),                                    # 9.8e-3
            Float16.(collect(-6.0:0.1:6.0)),                                    # 3.9e-2
            Float16.(@. 4sinh(0.5σ)/sinh(0.5)),                                 # 0.11, cells 12% apart
            (@. L*(s + 0.08sin(2π*s)/(2π))),                                    # 0.15
            (@. L*s + 1e-6*(L/32)*sin(2π*s)),                                   # 3.9e-7
        ]
        for z in uniform
            g = grid(z)
            @test g isa UniformGrid{float(eltype(z))}
            @test spacing(g) === float(eltype(z))((z[end] - z[1])/(length(z) - 1))
            @test nodes(g) == collect(float(eltype(z)), z)
            @test UniformGrid(z) isa UniformGrid
        end
        for z in nonuniform
            g = grid(z)
            @test g isa NonUniformGrid{float(eltype(z))}
            @test g.widths == cell_widths(nodes(g)) && narrowest(g) == minimum(g.widths)
            @test_throws ArgumentError UniformGrid(z)
        end
        # The refusal names the first node off.
        message = try
            UniformGrid([0.0, 1.0, 2.0, 3.5, 4.5])
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("at node 3", message)
        # Nodes are finite, increasing, at least two. Cell widths, which
        # `OnGrid` took until 0.2, are not nodes, and the refusal says so.
        @test_throws ArgumentError grid([1.0])
        @test_throws ArgumentError grid([0.0, NaN, 1.0])
        @test_throws ArgumentError grid([0.0, 1.0, Inf])
        @test_throws ArgumentError grid([0.0, 2.0, 1.0])
        message = try
            grid(fill(0.5, 8))
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("widths", message)
        # Integers are differenced in Float64: Int8 wrapped the span of
        # [-100, 0, 100] to -56, and the mean spacing with it.
        @test spacing(grid(Int8[-100, 0, 100])) === 100.0
        @test spacing(UniformGrid(Int16[-20000, 0, 20000])) === 20000.0
        message = try
            UniformGrid(Int8[-100, -99, 100])
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("spacing 199.0 at node 2", message)
    end

    @testset "a vector's spacing is its mean, and its cells span its box" begin
        # n cells of the first difference carry the rounding of the largest
        # node: 2.3e-10 off the box of the grid near 1e6, the mean 7.4e-12, and
        # the widths sum to n·h to rounding. Measured on Julia 1.13.
        z = collect(1e6 .+ (0:63)*0.1)
        box = 64*0.1
        @test abs(64*(z[2] - z[1]) - box)/box > 1e-10
        @test abs(64*spacing(grid(z)) - box)/box < 1e-11
        @test sum(widths(grid(z))) ≈ 64*spacing(grid(z)) rtol = 1e-14
    end

    @testset "constructors and accessors" begin
        g = UniformGrid(0.1, 8)
        @test !hasnodes(g) && length(g) == 8 && narrowest(g) === 0.1
        @test_throws ArgumentError nodes(g)
        @test widths(g) == fill(0.1, 8) && cell_widths(g) == widths(g)
        @test eltype(UniformGrid(0.1f0, 4)) === Float32 && eltype(UniformGrid(1//5, 4)) === Float64
        for h in (0.0, -0.1, Inf, NaN)
            @test_throws ArgumentError UniformGrid(h, 8)
        end
        @test_throws ArgumentError UniformGrid(0.1, 0)
        # The spacing as given, beside nodes whose differences round otherwise,
        # if it is theirs: the field solve and the Courant numbers take it, the
        # sums the nodes.
        xc = [(j - 0.5)*4π/6 for j in 1:6]
        gc = UniformGrid(4π/6, xc)
        @test spacing(gc) === 4π/6 && nodes(gc) == xc && xc[2] - xc[1] != 4π/6
        @test_throws ArgumentError UniformGrid(0.1, [0.0, 0.2, 0.1])
        @test_throws ArgumentError UniformGrid(0.2, collect(0.1:0.1:1.0))
        @test_throws ArgumentError UniformGrid(4π/6*(1 + 1e-6), xc)
        @test_throws ArgumentError UniformGrid(0.1, [0.0, 0.1, 0.3])
        # Built through the parametric constructors too, a grid is checked in
        # its own type: a spacing or a width that rounds to 0 in it, nodes that
        # do not increase.
        @test_throws ArgumentError UniformGrid{Float32, Nothing}(1e-50, 8, nothing)
        @test_throws ArgumentError UniformGrid{Float64, Vector{Float64}}(0.1, 3, [0.0, 0.2, 0.1])
        @test_throws ArgumentError NonUniformGrid{Float32, Nothing}([1e-50, 1.0], 1e-50, nothing)
        @test_throws ArgumentError NonUniformGrid{Float64, Vector{Float64}}([1.0, 1.0], 1.0, [1.0, 1.0])
        # A non-uniform grid has widths, the narrowest, and no spacing.
        w = [0.5 + 0.3sin(2π*i/16) for i in 1:16]
        n = NonUniformGrid(; widths = w)
        @test widths(n) == w && narrowest(n) == minimum(w) && !hasnodes(n) && length(n) == 16
        @test_throws MethodError spacing(n)
        @test_throws ArgumentError NonUniformGrid(; widths = [0.1, 0.0, 0.1])
        @test_throws ArgumentError NonUniformGrid(; widths = Float64[])
        @test_throws ArgumentError NonUniformGrid{Float64, Nothing}([0.1, 0.2], 0.2, nothing)
        z = cumsum(w)
        @test NonUniformGrid(z).widths == cell_widths(z) && nodes(NonUniformGrid(z)) == z
        # `cell_widths` moved here, and the driver's module still exports it.
        @test Vasilek.VlasovPoisson1D1V.cell_widths === cell_widths
        @test cell_widths(Float32[0, 1, 3, 4]) == Float32[1, 1.5, 1.5, 1]
        @test occursin("h = 0.1", sprint(show, g)) && occursin("NonUniformGrid", sprint(show, n))
    end

    @testset "the same cells, to the rounding either grid carries" begin
        g = UniformGrid(0.1, 8)
        @test samecells(g, UniformGrid(0.1, 8)) && !samecells(g, UniformGrid(0.1, 9))
        @test samecells(g, UniformGrid(0.1*(1 + 1e-11), 8))      # 1e-10 of a spacing
        @test !samecells(g, UniformGrid(0.1*(1 + 1e-9), 8))
        # 4 ulps of the coarser type
        @test samecells(UniformGrid(0.1f0, 8), g)
        @test !samecells(UniformGrid(0.1f0*(1 + 1f-5), 8), g)
        # or the rounding of the nodes a spacing is a difference of: 2.3e-10 of
        # 0.1 near 1e6, its first difference
        far = grid(collect(1e6 .+ (0:63)*0.1))
        @test samecells(far, UniformGrid(0.1, 64)) && samecells(UniformGrid(0.1, 64), far)
        @test !samecells(far, UniformGrid(0.1*(1 + 1e-7), 64))
        # each in its own type: Float32 nodes from -L/2 on L = 20π, whose
        # first difference is 1.0e-6 of L/N off, are a grid of L/N; their
        # rounding taken in Float64 was lost, and an OnGrid on L/N refused them.
        let L = 20π, N = 64
            x₃₂ = Float32.(collect(range(-L/2; step = L/N, length = N)))
            @test samecells(grid(x₃₂), UniformGrid(L/N, N))
            og = OnGrid(Upwind(), L/N, N)
            @test OnGrid(og, x₃₂) === og
        end
        w = [0.5 + 0.3sin(2π*i/16) for i in 1:16]
        @test samecells(NonUniformGrid(; widths = w), NonUniformGrid(; widths = w .* (1 + 1e-13)))
        @test !samecells(NonUniformGrid(; widths = w), NonUniformGrid(; widths = reverse(w)))
        # A uniform grid and a non-uniform one never have, equal widths or not.
        @test !samecells(UniformGrid(0.5, 4), NonUniformGrid(; widths = fill(0.5, 4)))
    end

    @testset "OnGrid dispatches on the grid and on what the scheme takes" begin
        @test takes(Upwind()) === Courant() && takes(PFC(fmin = 0, fmax = 1)) === Courant()
        @test takes(PFCNonUniform(fill(0.1, 8); fmin = 0, fmax = 1)) === Displacement()
        # An OnGrid takes a displacement, whatever its scheme takes.
        @test takes(OnGrid(Upwind(), 0.1, 8)) === Displacement()
        # A scheme of one's own that takes a displacement goes on any grid, and
        # an OnGrid of it is held to its grid as one of a Courant scheme is.
        og = OnGrid(DisplacedShift(), UniformGrid(0.1, 8))
        @test OnGrid(og, UniformGrid(0.1, 8)) === og && substep_bound(og) === 0.1
        @test_throws ArgumentError OnGrid(og, UniformGrid(0.2, 8))
        @test_throws ArgumentError OnGrid(og, UniformGrid(0.1, 9))
        # A PFCNonUniform is held to its widths in the coarser type's 4 ulps,
        # not to the nodes' rounding, which admitted widths a 1e-10 apart and a
        # run that summed on other widths than it advected on.
        far = collect(1e6 .+ (0:63)*0.1)
        p = PFCNonUniform(cell_widths(far) .* (1 + 1e-10); fmin = 0.0, fmax = 1.0)
        @test_throws ArgumentError OnGrid(p, far)
        # Its nodes' cell widths are 9.2e-10 of a cell from the grid's spacing,
        # which the run sums with: refused too. Built on the grid, it is taken.
        @test_throws ArgumentError OnGrid(PFCNonUniform(cell_widths(far); fmin = 0.0, fmax = 1.0), far)
        @test OnGrid(PFCNonUniform(grid(far); fmin = 0.0, fmax = 1.0), far) isa OnGrid
        n = 64
        src = [1.0 + 0.5sin(2π*i/n) for i in 1:n]
        x = collect(range(0.1; step = 0.1, length = n))
        stretched = [0.1 + 0.03sin(2π*i/n) for i in 1:n]
        # Every pair a run puts together, inferred, and its step allocating
        # nothing, sub-steps included: three at α = 0.25.
        pairs = [(Upwind(), UniformGrid(0.1, n)),
                 (Upwind(), grid(x)),
                 (Upwind(), UniformGrid(0.1f0, n)),
                 (Godunov(PiecewiseLinear(), VanLeer()), UniformGrid(0.1, n)),
                 (SemiLagrangian(CubicSpline()), grid(x)),
                 (PFC(fmin = 0.0, fmax = 2.0), UniformGrid(0.1, n)),
                 (PFCNonUniform(fill(0.1, n); fmin = 0.0, fmax = 2.0), UniformGrid(0.1, n)),
                 (PFCNonUniform(cell_widths(x); fmin = 0.0, fmax = 2.0), grid(x)),
                 (PFCNonUniform(stretched; fmin = 0.0, fmax = 2.0), NonUniformGrid(; widths = stretched))]
        for (scheme, g) in pairs
            og = @inferred OnGrid(scheme, g)
            @test og.grid === g
            ws = workspace(og, n, Float64)
            dest = similar(src)
            for α in (0.05, 0.25)
                @inferred advect!(dest, src, og, α, ws)
                @test (@allocated advect!(dest, src, og, α, ws)) == 0
            end
        end
        # A scheme that takes a Courant number goes on a uniform grid only.
        for s in (Upwind(), PFC(fmin = 0, fmax = 1), SemiLagrangian(CubicSpline()))
            @test_throws ArgumentError OnGrid(s, NonUniformGrid(; widths = stretched))
        end
        # The longest single step: the spacing, for a scheme that takes a
        # Courant number; a PFCNonUniform's own narrowest width.
        @test substep_bound(OnGrid(Upwind(), 0.1, n)) === 0.1
        p = PFCNonUniform(stretched; fmin = 0.0, fmax = 2.0)
        @test substep_bound(OnGrid(p)) === p.Δxmin === minimum(stretched)
    end

    @testset "a PFCNonUniform on a uniform grid steps by its own narrowest width" begin
        # `cell_widths` of a collected uniform range undercut its first spacing
        # by an ulp; a sub-step of the grid's spacing would cross the scheme's
        # narrowest cell, which its `advect!` refuses. `OnGrid` bounds the
        # sub-steps by the scheme's `Δxmin`, as it did before the grids.
        x = collect(range(4π/32; step = 4π/32, length = 32))
        h = x[2] - x[1]
        p = PFCNonUniform(cell_widths(x); fmin = 0.0, fmax = 1.0)
        @test minimum(p.Δx) < h
        src = [0.5 + 0.4sin(2π*i/32) for i in 1:32]
        @test_throws DomainError advect!(similar(src), src, p, h, workspace(p, 32))
        og = OnGrid(p, x)
        @test og.grid isa UniformGrid && substep_bound(og) === p.Δxmin < spacing(og.grid)
        @test nsubsteps(og, 2h) == 3
        @test all(isfinite, advect!(similar(src), src, og, 2h, workspace(og, 32)))
    end

    @testset "a PFCNonUniform on a grid's widths" begin
        # On a uniform grid every width is the spacing, and the limiter
        # coefficient 2 to the bit.
        p = @inferred PFCNonUniform(UniformGrid(0.1, 8); fmin = 0.0, fmax = 1.0)
        @test p isa PFCNonUniform{Float64, true} && p.Δx == fill(0.1, 8) && all(==(2), p.ξ)
        w = [0.5 + 0.3sin(2π*i/16) for i in 1:16]
        @test PFCNonUniform(NonUniformGrid(; widths = w); fmin = 0, fmax = 1).Δx == w
        @test PFCNonUniform(UniformGrid(0.1f0, 8); fmin = 0, fmax = 1) isa PFCNonUniform{Float32}
        # One on the nodes' own widths is held to the grid's to 1e-12 in
        # Float64. They are an ulp of the largest node off: on [0, 4π) 1.8e-13 of
        # a cell at 4096 nodes, taken, and 1.6e-12 at 16384, refused with how far
        # and what to build instead.
        for (N, taken) in ((4096, true), (16384, false))
            x = collect((0:N-1) .* (4π/N))
            p = PFCNonUniform(cell_widths(x); fmin = 0, fmax = 1)
            if taken
                @test OnGrid(p, x).grid isa UniformGrid
            else
                message = try
                    OnGrid(p, x)
                    ""
                catch err
                    sprint(showerror, err)
                end
                @test occursin("of a cell off", message) &&
                      occursin("PFCNonUniform(Vasilek.Grids.grid(x)", message)
            end
        end
    end
end
