using Vasilek
using NumericalIntegration

"""
Free streaming: the one exact kinetic solution available without a field solve.

With `E = 0` the Vlasov equation is `∂f/∂t + v ∂f/∂x = 0`, whose solution is
`f₀(x − vt, v)` exactly. Starting from a Maxwellian modulated in space,

    f₀(x, v) = (2π)^(-1/2) exp(−v²/2) · (1 + α cos kx)

the density integrates to a mode whose amplitude decays as a **Gaussian in
time**:

    n(x, t) = 1 + α exp(−k²t²/2) cos kx

because `∫ (2π)^(-1/2) exp(−v²/2) cos(kvt) dv = exp(−k²t²/2)`. This is phase
mixing: nothing is dissipated, the mode is carried into ever finer structure in
`v` by the shear, and the velocity integral of that structure cancels.

On a **discrete** velocity grid the cancellation is not permanent. The `Nv`
sheared modes rephase, and the density mode returns essentially in full at

    T_R = 2π/(kΔv)

which is a property of the grid rather than of the physics, and a number no
scheme can fudge -- it is set by `Δv` alone.

**What this covers that nothing else does.** The x-sweep and the velocity moment
are exercised together, against an analytic answer, outside the extended gate --
today that path runs only inside `verification_harness.jl`, behind
`VASILEK_EXTENDED=1`. And the ranking it produces is an accuracy comparison on a
*physical* observable rather than on an L² norm against a shifted sine, which is
what the rest of the advection suite measures.

The Courant limit binds here in a way it does not elsewhere: the fastest row
carries `max|v|·Δt/Δx`, so the velocity window and the time step are not
independent. Everything below runs at 0.611.
"""

const FS_L  = 4π
const FS_K  = 2π/FS_L                      # = 0.5
const FS_NX = 64
const FS_ΔT = 0.02
const FS_α  = 1e-2

"Amplitude of the `cos kx` component of `n`, on the uniform grid `x`."
mode_amplitude(n, x, k) = 2*sum(n .* cos.(k .* x))/length(x)

"""
    free_stream(scheme, v, nsteps; fmax)

Advect a modulated Maxwellian in `x` only, returning the time series of the
density mode amplitude. One row per velocity, one shared scheme value and one
workspace -- the reentrant pattern `test_threading` pins.
"""
function free_stream(scheme, v, nsteps)
    Δx = FS_L/FS_NX
    x = [(j-1)*Δx for j = 1:FS_NX]
    fM = @. exp(-0.5*v^2)/sqrt(2π)
    f = [fM[i]*(1 + FS_α*cos(FS_K*x[j])) for i in eachindex(v), j = 1:FS_NX]

    buf = Vector{Float64}(undef, FS_NX)
    ws = workspace(scheme, FS_NX)
    amps = Vector{Float64}(undef, nsteps + 1)
    for s = 0:nsteps
        n = [integrate(v, @view f[:, j]) for j = 1:FS_NX]
        amps[s+1] = mode_amplitude(n, x, FS_K)
        s == nsteps && break
        for i in eachindex(v)
            row = @view f[i, :]
            advect!(buf, row, scheme, v[i]*FS_ΔT/Δx, ws)
            copyto!(row, buf)
        end
    end
    return amps
end

@testset "Free streaming" begin

    @testset "the density mode decays as exp(−k²t²/2)" begin
        # Measured relative error in the mode amplitude over t ≤ 4, on
        # v ∈ [−6, 6] with Δv = 0.1 (Courant 0.611 on the fastest row):
        #
        #   SemiLagrangian cubic  2.7e-6
        #   PFC                   3.3e-5
        #   LaxWendroff           1.2e-3
        #   Upwind                1.7e-2
        #
        # Four orders of magnitude between the ends, on a quantity with physical
        # meaning. The tolerances below are each about 3x the measured value.
        #
        # Comparing pointwise against the Gaussian at every step is what makes
        # this sharp: a scheme that damped the mode exponentially -- which is
        # what excessive numerical diffusion looks like -- tracks the analytic
        # curve at one time and misses it everywhere else.
        v = collect(-6:0.1:6)
        nsteps = round(Int, 4.0/FS_ΔT)
        t = [s*FS_ΔT for s = 0:nsteps]
        exact = @. FS_α*exp(-0.5*FS_K^2*t^2)

        cases = [("SemiLagrangian cubic", SemiLagrangian(CubicSpline()), 1e-5),
                 ("PFC",                  PFC(fmin = 0.0, fmax = 0.5),   1e-4),
                 ("LaxWendroff",          LaxWendroff(),                 4e-3),
                 ("Upwind",               Upwind(),                      6e-2)]

        errors = Dict{String,Float64}()
        for (name, scheme, tol) in cases
            amps = free_stream(scheme, v, nsteps)
            err = maximum(abs, amps .- exact)/FS_α
            errors[name] = err
            println("  ", rpad(name, 22), "max rel error over t ≤ 4 = ", err)
            @test err < tol
        end

        # The ordering itself, which is the comparison this file exists to make.
        # It is deterministic -- these are errors, not timings -- so it can be
        # asserted rather than reported.
        @test errors["SemiLagrangian cubic"] < errors["PFC"] <
              errors["LaxWendroff"] < errors["Upwind"]
    end

    @testset "the mode recurs at T_R = 2π/(kΔv)" begin
        # Deliberately on a *coarse* velocity grid. Recurrence is a property of
        # Δv, so coarsening it brings T_R forward -- 31.4 here against 125.7 at
        # Δv = 0.1 -- and makes the test four times cheaper without weakening
        # it. The decay test above is where a fine grid matters.
        #
        # Measured: both schemes peak at t = 31.42 against T_R = 31.4159, one
        # time step away, recovering 1.0000 (semi-Lagrangian) and 0.9991 (PFC)
        # of the initial amplitude. The recurrence is a grid artefact and no
        # scheme can move it; what a scheme *can* do is damp the amplitude it
        # returns with, which is why both are asserted.
        Δv = 0.4
        v = collect(-6:Δv:6)
        T_R = 2π/(FS_K*Δv)
        nsteps = round(Int, 1.3*T_R/FS_ΔT)

        for (name, scheme) in [("SemiLagrangian cubic", SemiLagrangian(CubicSpline())),
                               ("PFC", PFC(fmin = 0.0, fmax = 0.5))]
            amps = free_stream(scheme, v, nsteps)
            # search past the initial decay so the t = 0 peak is not found
            from = round(Int, 0.5*T_R/FS_ΔT)
            i = argmax(abs.(amps[from:end])) + from - 1
            t_peak = (i - 1)*FS_ΔT
            recovered = abs(amps[i])/FS_α
            println("  ", rpad(name, 22), "recurrence at t = ", round(t_peak; digits = 3),
                    " (T_R = ", round(T_R; digits = 4), "), recovered ",
                    round(recovered; digits = 4), " of the initial amplitude")
            @test isapprox(t_peak, T_R; atol = 2*FS_ΔT)
            @test recovered > 0.9
        end
    end
end
