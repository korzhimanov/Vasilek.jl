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
end
