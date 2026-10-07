using Vasilek
using Vasilek: StrangSplitting
using Vasilek.StrangSplitting: strang_step!, Collide
using Vasilek.Advection: nsubsteps
using FFTW
using LinearAlgebra: mul!

# The splitting is tested here *in isolation from the advection schemes*, by
# handing it an exact spectral shift as its advection operator. A scheme's
# spatial error does not vanish as Δt → 0 -- semi-Lagrangian interpolation error
# in fact accumulates with the step count -- so measuring the temporal order
# against a real scheme measures the scheme, not the splitting. With an exact
# shift the only error left is the splitting error, and it comes out at the
# second order Strang promises.
#
# Every step goes through `strang_step!` with the schemes on their grids,
# `OnGrid(scheme, Δz)`, as `vlasov_poisson` takes it: `f[x, v]`, `cx` the
# displacements of the x lines, `cv` handed the state transposed.

# Rigid-rotation geometry, shared by the two testsets at the foot of this file.
const LD_ROT = 12.0
const Σ_ROT = 0.8
const T_ROT = 1.0

"""
    SpectralShift()

An exact periodic translation, through the spectrum: `advect!` shifts a line
by `c` cells, `c` being its fourth argument -- the Courant number an `OnGrid`
divides a displacement into. No Courant limit: like `SemiLagrangian`, it takes
a step of any width whole.
"""
struct SpectralShift <: AbstractAdvection1D end

struct SpectralShiftWorkspace{P, Q}
    in::Vector{Float64}
    out::Vector{Float64}
    spectrum::Vector{ComplexF64}
    k::Vector{Float64}          # radians per cell
    forward::P
    inverse::Q
end

function Vasilek.workspace(::SpectralShift, n::Integer, ::Type = Float64)
    in = Vector{Float64}(undef, n)
    spectrum = Vector{ComplexF64}(undef, n ÷ 2 + 1)
    return SpectralShiftWorkspace(in, similar(in), spectrum, collect(2π .* rfftfreq(n)),
                                  plan_rfft(in), plan_irfft(copy(spectrum), n))
end

function Vasilek.Advection.advect!(dest, src, ::SpectralShift, c, ws::SpectralShiftWorkspace)
    # Through the workspace's own arrays: a plan is made for an aligned,
    # contiguous one and refuses a column view that is not aligned as it was, as
    # `PoissonFourier1D.solve!` copies for the same reason.
    copyto!(ws.in, src)
    mul!(ws.spectrum, ws.forward, ws.in)
    @. ws.spectrum *= cis(-ws.k*c)
    mul!(ws.out, ws.inverse, ws.spectrum)
    return copyto!(dest, ws.out)
end

Vasilek.Advection.nsubsteps(::OnGrid{SpectralShift}, α) = 1

"""
A collision operator that changes nothing and records the `Δt` of every call:
which lines the hook collides, how often, and for how long.
"""
struct CountingCollisions <: AbstractCollisionOperator
    Δts::Vector{Float64}
end

function Vasilek.Collisions.collide!(dest, src, op::CountingCollisions, v, Δt, ::Nothing)
    push!(op.Δts, Δt)
    return copyto!(dest, src)
end

"""
    written_out!(g, sx, sv, cx, cv[, op, v, Δt])

One Strang step of `g[x, v]` written out from `advect!` and `copyto!` as the 0.1
form took it -- the first direction over the columns of `g`, the second over
those of its transpose, `cv` handed that transpose -- and with `op`, every line
of the second direction collided for `Δt/2` either side of its advection.
"""
function written_out!(g, sx, sv, cx, cv, op = nothing, v = nothing, Δt = nothing)
    nx, nv = size(g)
    bx, bv = zeros(nx), zeros(nv)
    wx, wv = workspace(sx, nx), workspace(sv, nv)
    wc = op === nothing ? nothing : workspace(op, nv, Float64)
    for j in 1:nv
        advect!(bx, view(g, :, j), sx, cx[j]/2, wx)
        g[:, j] = bx
    end
    f = Matrix(g')
    c = cv(f)
    for i in 1:nx
        line = view(f, :, i)
        op === nothing || (collide!(bv, line, op, v, Δt/2, wc); copyto!(line, bv))
        advect!(bv, line, sv, c[i], wv)
        copyto!(line, bv)
        op === nothing || (collide!(bv, line, op, v, Δt/2, wc); copyto!(line, bv))
    end
    g .= f'
    for j in 1:nv
        advect!(bx, view(g, :, j), sx, cx[j]/2, wx)
        g[:, j] = bx
    end
    return g
end

@testset "StrangSplitting" begin

    @testset "a zero step is the identity" begin
        # With no displacement in either direction the three sweeps and the two
        # transposes must leave the data alone. Deliberately non-square: Nx ≠ Nv
        # is the case the transposes can get wrong.
        Nx, Nv = 32, 24
        f₀ = [1.0 + 0.3*sin(2π*i/Nx)*cos(2π*j/Nv) for i = 1:Nx, j = 1:Nv]
        for scheme in (SemiLagrangian(CubicSpline()), SpectralShift())
            sx, sv = OnGrid(scheme, fill(1.0, Nx)), OnGrid(scheme, fill(1.0, Nv))
            f = copy(f₀)
            @test strang_step!(f, sx, sv, zeros(Nv), _ -> zeros(Nx), workspace(sx, sv, f)) === f
            println("  zero step, ", rpad(nameof(typeof(scheme)), 15),
                    "max|f - f₀| = ", maximum(abs, f .- f₀))
            @test size(f) == (Nx, Nv)
            @test maximum(abs, f .- f₀) ≤ 1e-13
        end
    end

    @testset "second order in Δt" begin
        # Rigid rotation in phase space: ∂f/∂t + v ∂f/∂x − x ∂f/∂v = 0, whose
        # exact solution is the initial condition rotated by −t. Both the
        # split-step reference and the analytic answer are checked, because
        # agreeing with a fine reference only proves self-consistency.
        #
        # Measured orders: 2.003, 2.001, 2.002 against the analytic rotation.
        N = 64; L = 12.0; h = L/N; σ = 0.8; T = 1.0
        z = [-L/2 + (i-1)*h for i = 1:N]        # x and v share this grid
        s = OnGrid(SpectralShift(), fill(h, N))

        blob(a, b) = [exp(-((z[i]-a)^2 + (z[j]-b)^2)/(2σ^2)) for i = 1:N, j = 1:N]

        function rotate(Δt)
            f = blob(2.0, 0.0)
            ws = workspace(s, s, f)
            for _ = 1:round(Int, T/Δt)
                strang_step!(f, s, s, z .* Δt, _ -> -z .* Δt, ws)
            end
            return f
        end

        exact = blob(2.0*cos(T), -2.0*sin(T))
        errs = Float64[]
        for Δt in (T/8, T/16, T/32, T/64)
            e = maximum(abs, rotate(Δt) .- exact)
            isempty(errs) || println("  Δt = ", rpad(round(Δt; digits = 6), 10),
                                     " err = ", rpad(round(e; sigdigits = 4), 11),
                                     " order = ", round(log2(errs[end]/e); digits = 3))
            push!(errs, e)
        end
        for i = 2:length(errs)
            @test isapprox(log2(errs[i-1]/errs[i]), 2.0; atol = 0.15)
        end
        # and the absolute error at the coarsest step, so a uniform inflation
        # that preserves the slope still fails
        @test errs[1] < 5e-3
    end

    @testset "strang_step! is the composition, bit for bit" begin
        # `strang_step!` has to take the step `written_out!` takes, to the bit:
        # its transposes, its half steps of `cx[j]*(1//2)` against `cx[j]/2`,
        # and `cv` handed the transposed state. For a uniform scheme, for
        # PFCNonUniform on grids that are not uniform, and with displacements
        # wide enough that `OnGrid` splits them.
        Nx, Nv = 32, 24
        g₀ = [exp(-((i-16)/5)^2 - ((j-12)/4)^2) for i = 1:Nx, j = 1:Nv]
        cx = collect(range(-1.6, 1.6, length = Nv))
        cv(ft) = [0.4*sin(2π*i/Nx) + 1e-3*sum(view(ft, :, i)) for i in 1:Nx]
        wx = [0.5 + 0.2sin(2π*i/Nx) for i in 1:Nx]
        wv = [0.25 + 0.1cos(2π*j/Nv) for j in 1:Nv]
        for (sx, sv) in ((OnGrid(SemiLagrangian(CubicSpline()), fill(0.5, Nx)),
                          OnGrid(SemiLagrangian(CubicSpline()), fill(0.25, Nv))),
                         (OnGrid(PFC(fmin = 0.0, fmax = 1.0), fill(0.5, Nx)),
                          OnGrid(Godunov(PiecewiseLinear(), VanLeer()), fill(0.25, Nv))),
                         (OnGrid(PFCNonUniform(wx; fmin = 0.0, fmax = 1.0)),
                          OnGrid(PFCNonUniform(wv; fmin = 0.0, fmax = 1.0))))
            g = copy(g₀)
            for _ in 1:5
                written_out!(g, sx, sv, cx, cv)
            end
            f = copy(g₀)
            ws = workspace(sx, sv, f)
            for _ in 1:5
                @test strang_step!(f, sx, sv, cx, cv, ws) === f
            end
            @test f == g
            sx.scheme isa SemiLagrangian ||
                @test nsubsteps(sx, cx[end]/2) > 1 && nsubsteps(sv, 0.4) > 1
        end
    end

    @testset "OnGrid splits a step wider than a cell into equal ones" begin
        # `m` sub-steps of `α/m`, each taken as the scheme takes it, to the bit;
        # within a cell, the scheme's own single step. `h` is the narrowest
        # cell, 0.2 on the stretched grid.
        n = 40
        src = [exp(-((i - 15)/4)^2) for i in 1:n]
        uniform = fill(0.5, n)
        stretched = [0.5 + 0.3sin(2π*i/n) for i in 1:n]
        for (scheme, Δz) in ((Upwind(), uniform), (Godunov(PiecewiseLinear(), VanLeer()), uniform),
                             (PFCNonUniform(stretched; fmin = 0.0, fmax = 1.0), stretched))
            og = OnGrid(scheme, Δz)
            ws, inner = workspace(og, n), workspace(scheme, n)
            h = minimum(Δz)
            argument(a) = scheme isa PFCNonUniform ? a : a/h
            for α in (2.5h, -2.5h, 0.7h, -h)
                m = nsubsteps(og, α)
                @test m == (abs(α) > h ? 3 : 1)
                ref, buf = copy(src), similar(src)
                for _ in 1:m
                    advect!(buf, ref, scheme, argument(α/m), inner)
                    copyto!(ref, buf)
                end
                @test advect!(similar(src), src, og, α, ws) == ref
            end
            # the four-argument form builds its own workspace
            @test advect!(similar(src), src, og, 2.5h) == advect!(similar(src), src, og, 2.5h, ws)
        end
        # SemiLagrangian has no Courant limit, and takes the step whole.
        s = SemiLagrangian(CubicSpline())
        og = OnGrid(s, uniform)
        @test nsubsteps(og, 1.25) == 1
        @test advect!(similar(src), src, og, 1.25, workspace(og, n)) ==
              advect!(similar(src), src, s, 2.5, workspace(s, n))
    end

    @testset "OnGrid refuses a grid its scheme cannot take" begin
        # A uniform scheme on a non-uniform grid would take one of its spacings
        # and be wrong by the ratio between them.
        @test_throws ArgumentError OnGrid(Upwind(), [1.0, 2.0, 1.0, 1.0])
        message = try
            OnGrid(LaxWendroff(), [1.0, 2.0, 1.0, 1.0])
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("LaxWendroff", message) && occursin("PFCNonUniform", message)
        @test OnGrid(Upwind(), [1.0, 1.0 + 1e-14, 1.0]).h == 1.0     # uniform to 1e-12
        # A PFCNonUniform carries its grid, and is refused on another.
        Δz = [0.5 + 0.3sin(2π*i/16) for i in 1:16]
        p = PFCNonUniform(Δz; fmin = 0.0, fmax = 1.0)
        @test OnGrid(p).h == OnGrid(p, Δz).h == minimum(Δz)
        @test_throws ArgumentError OnGrid(p, reverse(Δz))
        @test_throws DimensionMismatch OnGrid(p, Δz[1:end-1])
        @test_throws ArgumentError OnGrid(p, 0.5, 16)                # not its narrowest cell
        # Already on a grid: kept on that grid, refused on another.
        og = OnGrid(Upwind(), fill(0.5, 8))
        @test OnGrid(og, fill(0.5, 8)) === og
        @test_throws ArgumentError OnGrid(og, fill(0.25, 8))
        @test_throws ArgumentError OnGrid(og, 0.5, 8)
        # To the same 1e-12 the grid is taken at: `L/N` and the widths of the
        # cell centres `(j - 1/2)L/N` differ in the last bit (L = 4π, N = 6), and
        # an exact comparison refused an OnGrid on the very grid it was given.
        L, N = 4π, 6
        widths = Vasilek.VlasovPoisson1D1V.cell_widths([(j - 0.5)*L/N for j in 1:N])
        og₆ = OnGrid(Upwind(), fill(L/N, N))
        @test first(widths) != L/N
        @test OnGrid(og₆, widths) === og₆
        # No cells, or cells of no positive width.
        @test_throws ArgumentError OnGrid(Upwind(), Float64[])
        @test_throws ArgumentError OnGrid(Upwind(), fill(-0.5, 8))
        # The grid's length binds the workspace and the line.
        @test_throws DimensionMismatch workspace(og, 9)
        @test_throws DimensionMismatch advect!(zeros(9), ones(9), og, 0.25)
    end

    @testset "the hook collides half a step either side of the kick" begin
        # `X(Δt/2) · C(Δt/2) K(Δt) C(Δt/2) · X(Δt/2)`, against the composition
        # written out, bit for bit: two beams, so that BGK moves every line, on
        # the driver's default schemes.
        Nx, Nv = 16, 49
        v = collect(range(-6.0, 6.0; length = Nv))
        Δt = 0.2
        f₀ = [(exp(-(u - 1.5)^2/2) + exp(-(u + 1.5)^2/2))/(2sqrt(2π))*(1 + 0.3cos(2π*i/Nx))
              for i in 1:Nx, u in v]
        sx = OnGrid(PFCNonUniform(fill(0.5, Nx); fmin = 0.0, fmax = Inf))
        sv = OnGrid(PFCNonUniform(fill(v[2] - v[1], Nv); fmin = 0.0, fmax = Inf))
        cx = v .* Δt
        cv(ft) = [(0.3sin(2π*i/Nx) + 1e-3*sum(view(ft, :, i)))*Δt for i in 1:Nx]
        op = BGK(0.5)
        g = copy(f₀)
        for _ in 1:4
            written_out!(g, sx, sv, cx, cv, op, v, Δt)
        end
        f, plain = copy(f₀), copy(f₀)
        ws = workspace(sx, sv, f, op)
        for _ in 1:4
            strang_step!(f, sx, sv, cx, cv, ws, Collide(op, v, Δt))
            strang_step!(plain, sx, sv, cx, cv, ws)
        end
        @test f == g
        @test maximum(abs, f .- plain) > 1e-3*maximum(f)

        # Every line of the second direction, twice a step, for half of it each
        # time, and nothing else: an operator that changes nothing leaves the
        # step without collisions to the bit.
        counting = CountingCollisions(Float64[])
        f = copy(f₀)
        strang_step!(f, sx, sv, cx, cv, workspace(sx, sv, f, counting), Collide(counting, v, Δt))
        @test length(counting.Δts) == 2Nx
        @test all(==(Δt/2), counting.Δts)
        plain = copy(f₀)
        @test f == strang_step!(plain, sx, sv, cx, cv, workspace(sx, sv, plain))

        # A workspace built for no operator, or for another kind, is refused
        # before the first sweep, and leaves `f` as it was; it was a MethodError
        # from inside the kick, half a step in. The operator's parameters aside:
        # `BGK(0.5f0)`'s workspace serves `BGK(0.5)`.
        f = copy(f₀)
        message = try
            strang_step!(f, sx, sv, cx, cv, workspace(sx, sv, f), Collide(op, v, Δt))
            ""
        catch err
            sprint(showerror, err)
        end
        @test occursin("BGK", message) && occursin("no collision operator", message)
        @test_throws ArgumentError strang_step!(f, sx, sv, cx, cv, workspace(sx, sv, f, counting),
                                                Collide(op, v, Δt))
        @test f == f₀
        g = copy(f₀)
        strang_step!(f, sx, sv, cx, cv, workspace(sx, sv, f, BGK(0.5f0)), Collide(op, v, Δt))
        @test f == strang_step!(g, sx, sv, cx, cv, workspace(sx, sv, g, op), Collide(op, v, Δt))
    end

    @testset "strang_step! checks its shapes" begin
        h = rand(8, 6)
        s = Upwind()
        ws = workspace(s, s, h)
        @test_throws DimensionMismatch strang_step!(h, s, s, zeros(5), _ -> zeros(8), ws)
        @test_throws DimensionMismatch strang_step!(h, s, s, zeros(6), _ -> zeros(7), ws)
        @test_throws DimensionMismatch strang_step!(h, s, s, zeros(6), _ -> zeros(8),
                                                    workspace(s, s, rand(6, 8)))
        wc = workspace(s, s, h, BGK(1.0))
        @test_throws DimensionMismatch strang_step!(h, s, s, zeros(6), _ -> zeros(8), wc,
                                                    Collide(BGK(1.0), zeros(5), 0.1))
    end

    # ------------------------------------------------------------------
    # Everything above measures the splitting with the schemes deliberately
    # removed, by handing it an exact spectral shift. That is the right way to
    # measure the splitting *alone*, and it should stay -- but it means nothing
    # here has ever measured the splitting and a real scheme together against an
    # analytic answer, which is what a production run actually is.
    #
    # The two testsets below do that, on the same rigid rotation, so the
    # comparison against the exact-shift result above is direct.

    "The rotation grid: `z` is shared by x and v, and `h` is its spacing."
    rot_grid(N) = (LD_ROT/N, [-LD_ROT/2 + (i-1)*(LD_ROT/N) for i = 1:N])

    "A Gaussian blob at `(a, b)` in phase space."
    rot_blob(z, a, b) = [exp(-((z[i]-a)^2 + (z[j]-b)^2)/(2*Σ_ROT^2))
                         for i in eachindex(z), j in eachindex(z)]

    """
        rotate(scheme, N, nsteps; g0)

    `nsteps` Strang steps of `∂f/∂t + v ∂f/∂x − x ∂f/∂v = 0` over `T_ROT`, with
    `scheme` doing both directions.

    `nsteps = N` throughout, which is not tidiness: the rotation's displacement
    grows with `|z|`, so the fastest line runs at `(L/2)·Δt/h = N·T/(2·nsteps)`,
    and `PFC` is a finite-volume scheme with a Courant limit of 1. At
    `nsteps = N` that is 0.5; at the `nsteps = N/4` first tried it is 2.0, and
    `PFC` undershot to -3.9e-10 and tripped its own bounds check (`OnGrid` would
    now split those steps). Tying the step count to the resolution also means
    Δx, Δv and Δt refine together, so the orders below are joint rather than
    temporal.
    """
    function rotate(scheme, N, nsteps; g0 = nothing)
        h, z = rot_grid(N)
        Δt = T_ROT/nsteps
        f = g0 === nothing ? rot_blob(z, 2.0, 0.0) : copy(g0)
        s = OnGrid(scheme, fill(h, N))
        ws = workspace(s, s, f)
        for _ = 1:nsteps
            strang_step!(f, s, s, z .* Δt, _ -> -z .* Δt, ws)
        end
        return f
    end

    @testset "splitting and scheme together reproduce the analytic rotation" begin
        # Measured, refining Δx, Δv and Δt together over N = 32, 64, 128:
        #
        #   SemiLagrangian cubic   6.37e-3  5.24e-4  5.72e-5   orders 3.60 3.19
        #   PFC                    6.36e-2  1.33e-2  2.85e-3   orders 2.26 2.22
        #   LaxWendroff            1.17e-1  3.27e-2  8.39e-3   orders 1.84 1.96
        #
        # Each scheme keeps its own spatial order rather than being dragged to
        # the splitting's second: the cubic spline reaches 3.19 where Strang
        # alone would give 2. That is not an accident of this problem being
        # easy -- it is that a rigid rotation factors into shears, which is
        # exactly what the splitting does, so the commutator error is unusually
        # small here and the scheme is what is left. Worth knowing when reading
        # the second-order result the exact-shift test above reports: that is
        # the splitting's order, and this is what a real run sees.
        for (name, scheme, expected) in
                (("SemiLagrangian cubic", SemiLagrangian(CubicSpline()), 3.0),
                 ("PFC",                  PFC(fmin = -1e-12, fmax = 1.0), 2.0),
                 ("LaxWendroff",          LaxWendroff(),                  2.0))
            errs = Float64[]
            for N in (32, 64, 128)
                _, z = rot_grid(N)
                exact = rot_blob(z, 2*cos(T_ROT), -2*sin(T_ROT))
                push!(errs, maximum(abs, rotate(scheme, N, N) .- exact))
            end
            orders = [log2(errs[i-1]/errs[i]) for i in 2:length(errs)]
            println("  ", rpad(name, 22), join([string(round(e; sigdigits = 3), " ")
                                                for e in errs]),
                    "  orders ", join(round.(orders; digits = 2), " "),
                    "  (expected ", expected, ")")
            @test issorted(errs; rev = true)
            @test isapprox(orders[end], expected; atol = 0.3)
        end
    end

    @testset "the splitting is reversible, and dissipation is what breaks it" begin
        # Rotation is invariant under `(t, v) → (−t, −v)`, so rotating forward,
        # flipping the velocity axis, rotating forward again and flipping back
        # must return the initial state. Strang splitting is symmetric, so the
        # composition preserves this exactly; what does not is the scheme's own
        # dissipation, which has no time-reverse.
        #
        # That makes the round-trip error a direct measure of irreversibility,
        # and it separates the schemes far more sharply than the forward error
        # does. Measured at N = 64 and 128:
        #
        #   SemiLagrangian cubic   1.01e-3 -> 1.23e-4   (x8.2)
        #   LaxWendroff            3.77e-3 -> 4.94e-4   (x7.6)
        #   PFC                    1.95e-2 -> 3.61e-3   (x5.4)
        #   Upwind                 4.04e-1 -> 2.56e-1   (x1.6)
        #
        # Upwind's 0.40 is 40% of the peak: it has smeared the blob so far that
        # there is nothing left to reverse, and refining barely helps because
        # the dissipation is first order. It is included precisely because it
        # fails -- a test where every scheme passes would not show that this
        # measures dissipation rather than the splitting.
        #
        # The flip is `[1; N:-1:2]`, not `reverse`: on a periodic grid starting
        # at `-L/2`, index 1 is its own mirror (`-L/2 ≡ L/2`) and the rest
        # reverse around it. Plain `reverse` maps `z → -z - h` and is off by a
        # cell, which shows up as a first-order error that no refinement clears.
        flip(A, N) = A[:, [1; N:-1:2]]

        for (name, scheme, tol) in
                (("SemiLagrangian cubic", SemiLagrangian(CubicSpline()),  3e-4),
                 ("LaxWendroff",          LaxWendroff(),                  1e-3),
                 ("PFC",                  PFC(fmin = -1e-12, fmax = 1.0), 1e-2))
            errs = Float64[]
            for N in (64, 128)
                _, z = rot_grid(N)
                g0 = rot_blob(z, 2.0, 0.0)
                there = rotate(scheme, N, N; g0 = g0)
                back = flip(rotate(scheme, N, N; g0 = flip(there, N)), N)
                push!(errs, maximum(abs, back .- g0))
            end
            println("  ", rpad(name, 22), "round trip ", round(errs[1]; sigdigits = 3),
                    " -> ", round(errs[2]; sigdigits = 3),
                    "  (x", round(errs[1]/errs[2]; digits = 1), ")")
            @test errs[2] < tol
            @test errs[2] < errs[1]/2        # refinement genuinely recovers it
        end

        # And the counter-example, which is what gives the three above their
        # meaning. First-order dissipation is not recoverable by refining at
        # this rate: measured 0.404 falling only to 0.256.
        errs = Float64[]
        for N in (64, 128)
            _, z = rot_grid(N)
            g0 = rot_blob(z, 2.0, 0.0)
            there = rotate(Upwind(), N, N; g0 = g0)
            back = flip(rotate(Upwind(), N, N; g0 = flip(there, N)), N)
            push!(errs, maximum(abs, back .- g0))
        end
        println("  ", rpad("Upwind", 22), "round trip ", round(errs[1]; sigdigits = 3),
                " -> ", round(errs[2]; sigdigits = 3),
                "  (x", round(errs[1]/errs[2]; digits = 1), ", and still 26% of the peak)")
        @test errs[2] > 0.1
        @test errs[1]/errs[2] < 2.5
    end
end
