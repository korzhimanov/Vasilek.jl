using Vasilek
using Vasilek: StrangSplitting
using FFTW
using LinearAlgebra: mul!

# `StrangSplitting` had no test of its own.
#
# It is exercised only by `verification_harness.jl`, which runs behind
# `VASILEK_EXTENDED=1` -- so in a default test run the module was executed zero
# times, and the transpose bookkeeping that ties the two phase-space layouts
# together was never checked at all.
#
# The splitting is tested here *in isolation from the advection schemes*, by
# handing it an exact spectral shift as its advection operator. A scheme's
# spatial error does not vanish as Δt → 0 -- semi-Lagrangian interpolation error
# in fact accumulates with the step count -- so measuring the temporal order
# against a real scheme measures the scheme, not the splitting. With an exact
# shift the only error left is the splitting error, and it comes out at the
# second order Strang promises.

# Rigid-rotation geometry, shared by the two testsets at the foot of this file.
const LD_ROT = 12.0
const Σ_ROT = 0.8
const T_ROT = 1.0

"Exact periodic translation by a physical displacement `α`, via the spectrum."
function exact_shift!(col, α, plan, iplan, k, buf)
    mul!(buf, plan, col)
    @. buf *= cis(-k*α)
    mul!(col, iplan, buf)
    return col
end

"An advection closure of the shape `make_time_step_2d!` expects."
function shift_advector(n, h)
    k = 2π*collect(rfftfreq(n, 1/h))
    plan = plan_rfft(Vector{Float64}(undef, n))
    buf = Vector{ComplexF64}(undef, length(k))
    iplan = plan_irfft(copy(buf), n)
    return (col, α) -> exact_shift!(col, α, plan, iplan, k, buf)
end

"One scheme-based advector per direction, as the solvers actually use it."
function scheme_advector(scheme, n, h)
    ws = workspace(scheme, n)
    buf = Vector{Float64}(undef, n)
    return (col, α) -> (advect!(buf, col, scheme, α/h, ws); copyto!(col, buf))
end

@testset "StrangSplitting" begin

    @testset "a zero step is the identity" begin
        # With no displacement in either direction the three sweeps must leave
        # the data alone, and the two layouts must still be transposes.
        # Deliberately non-square: Nx ≠ Nv is the case the `f[2][:] = (f[1])'`
        # bookkeeping can get wrong, and the only caller that exercises it is
        # behind the extended gate.
        Nx, Nv = 32, 24
        g = [1.0 + 0.3*sin(2π*i/Nx)*cos(2π*j/Nv) for i = 1:Nx, j = 1:Nv]
        f = Matrix(g')
        g₀ = copy(g)

        s = SemiLagrangian(CubicSpline())
        ax! = scheme_advector(s, Nx, 1.0)
        av! = scheme_advector(s, Nv, 1.0)

        StrangSplitting.make_time_step_2d!((g, f), (_ -> zeros(Nv), _ -> zeros(Nx)),
                                           (ax!, av!))
        println("  zero step: max|g - g₀| = ", maximum(abs, g .- g₀),
                ", max|f - g'| = ", maximum(abs, f .- Matrix(g')))
        @test maximum(abs, g .- g₀) ≤ 1e-13
        @test size(g) == (Nx, Nv) && size(f) == (Nv, Nx)
        @test maximum(abs, f .- Matrix(g')) ≤ 1e-13
    end

    @testset "second order in Δt" begin
        # Rigid rotation in phase space: ∂f/∂t + v ∂f/∂x − x ∂f/∂v = 0, whose
        # exact solution is the initial condition rotated by −t. Both the
        # split-step reference and the analytic answer are checked, because
        # agreeing with a fine reference only proves self-consistency.
        #
        # Measured orders: 2.003, 2.001, 2.000, 2.001 against the analytic
        # rotation, and the same against a Δt/1024 reference.
        N = 64; L = 12.0; h = L/N; σ = 0.8; T = 1.0
        z = [-L/2 + (i-1)*h for i = 1:N]        # x and v share this grid
        adv! = shift_advector(N, h)

        blob(a, b) = [exp(-((z[i]-a)^2 + (z[j]-b)^2)/(2σ^2)) for i = 1:N, j = 1:N]

        function rotate(Δt)
            g = blob(2.0, 0.0); f = Matrix(g')
            for _ = 1:round(Int, T/Δt)
                StrangSplitting.make_time_step_2d!(
                    (g, f), (_ -> z.*Δt, _ -> -z.*Δt), (adv!, adv!))
            end
            return g
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

    @testset "f[2] is stale on return" begin
        # The third sweep writes f[1] and never transposes back, so on return
        # f[1] is current and f[2] lags by the final half-step. Every caller in
        # this repository reads diagnostics off f[2], so the discrepancy is
        # real and worth pinning rather than discovering later.
        Nx, Nv = 32, 24
        g = [exp(-((i-16)/5)^2 - ((j-12)/4)^2) for i = 1:Nx, j = 1:Nv]
        f = Matrix(g')
        s = SemiLagrangian(CubicSpline())
        StrangSplitting.make_time_step_2d!(
            (g, f), (_ -> fill(0.3, Nv), _ -> fill(0.2, Nx)),
            (scheme_advector(s, Nx, 1.0), scheme_advector(s, Nv, 1.0)))
        lag = maximum(abs, f .- Matrix(g'))
        println("  max|f[2] - f[1]'| after a real step = ", lag,
                "   (nonzero: f[2] lags by the last half-step)")
        @test lag > 1e-6
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
    `PFC` undershoots to -3.9e-10 and trips its own bounds check. Tying the step
    count to the resolution also means Δx, Δv and Δt refine together, so the
    orders below are joint rather than temporal.
    """
    function rotate(scheme, N, nsteps; g0 = nothing)
        h, z = rot_grid(N)
        Δt = T_ROT/nsteps
        g = g0 === nothing ? rot_blob(z, 2.0, 0.0) : copy(g0)
        f = Matrix(g')
        ax! = scheme_advector(scheme, N, h)
        av! = scheme_advector(scheme, N, h)
        for _ = 1:nsteps
            StrangSplitting.make_time_step_2d!((g, f), (_ -> z.*Δt, _ -> -z.*Δt),
                                               (ax!, av!))
        end
        return g
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
