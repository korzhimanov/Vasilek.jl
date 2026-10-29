@isdefined(VM_HARNESS) || include(joinpath(@__DIR__, "vm_harness.jl"))

using Vasilek.VlasovMaxwell1D1P: distribution, density, kinetic_energy, longitudinal_energy,
                                 energy_account, step!
const VMM = Vasilek.VlasovMaxwell1D1P

# The relativistic Vlasov–Maxwell solver: what holds by construction is held to
# round-off here, on a foil small enough to run in seconds; the accuracy claims
# of `docs/src/verification.md` are in the extended block at the end.

# A foil of `l_nm` on a grid of `dx_nm`, a few micrometres either side, on
# reduced momentum grids, for the contracts.
small_foil(; kw...) = laser_foil(; x_um = (-2.0, 2.0), plasma_um = (-0.5, 0.6), pml_um = 0.4,
                                 electron_grid = (pmin = -2.0, pmax = 2.0),
                                 ion_grid = (pmin = -0.004, pmax = 0.004), kw...)

charge_density(vm) = sum(s.q .* density(vm, i) for (i, s) in enumerate(vm.species))

@testset "VlasovMaxwell1D1P" begin
    @testset "momentum grid and Maxwellian" begin
        e = momentum_edges(; pmin = -3.0, pmax = 7.0, Δp = 0.01, core = 0.5, ratio = 1.05, Δpmax = 0.2)
        Δ = diff(e)
        @test e[1] ≤ -3 && e[end] ≥ 7
        @test any(iszero, e)                          # an edge at 0
        @test all(≈(0.01), Δ[abs.(e[1:end-1]) .< 0.5 .- 1e-12])
        @test maximum(Δ) ≤ 0.2 + 1e-12
        @test all(abs.(Δ[2:end] ./ Δ[1:end-1] .- 1) .≤ 0.05 + 1e-9)   # growth at most 5%
        θ = 1000/MEC2_EV
        f = maxwellian(e, θ; n = 3.0)
        @test sum(f .* Δ) ≈ 3.0 rtol = 1e-14
        pc = (e[1:end-1] .+ e[2:end])./2
        # the temperature of the cell averages, `⟨p²⟩ = θ(1 + O(θ))` plus the
        # cell's own Δp²/12
        @test sum(f .* pc.^2 .* Δ)/3 ≈ θ rtol = 0.03
    end

    # `step!` allocating nothing is held in `test_allocations.jl`.
    @testset "a field-free foil: charge, Gauss's law" begin
        vm, meta = small_foil(; a0 = 0.0)
        N0 = [sum(density(vm, i))*vm.Δx for i in eachindex(vm.species)]
        for _ in 1:402
            step!(vm)
        end
        for (i, s) in enumerate(vm.species)
            N = sum(density(vm, i))*vm.Δx + s.lost_x + s.lost_p + s.flushed
            @test abs(N - N0[i])/N0[i] < 1e-12
        end
        # Gauss's law between every pair of faces, to round-off
        ρ = charge_density(vm)
        gauss = [vm.ex[c+1] - vm.ex[c] - vm.Δx*ρ[c] for c in 1:vm.Nc]
        println("  Gauss's law: max residual ", maximum(abs, gauss), " against max|Eₓ| ",
                maximum(abs, vm.ex))
        @test maximum(abs, gauss) < 1e-11*max(1, maximum(abs, vm.ex))
        @test minimum(distribution(vm, 1)) ≥ 0
    end

    # Bit for bit, the tallies of what left the grid included: the default
    # number of chunks follows the threads, and so would any difference.
    @testset "chunks do not change the answer, $scheme" for scheme in (PFCScheme(), SplineScheme())
        tallies(vm) = [(s.lost_x, s.lost_x_energy, s.lost_p, s.lost_p_energy, s.flushed)
                       for s in vm.species]
        a = small_foil(; a0 = 2.0, ξ = 0.3, scheme, chunks = 1)[1]
        b = small_foil(; a0 = 2.0, ξ = 0.3, scheme, chunks = 7)[1]
        for _ in 1:300
            step!(a); step!(b)
        end
        @test distribution(a, 1) == distribution(b, 1)
        @test a.ex == b.ex
        @test a.mesh.ey == b.mesh.ey
        @test tallies(a) == tallies(b)
        @test any(t -> t[3] != 0 && t[5] != 0, tallies(a))    # there was something to add
    end

    @testset "the spline's momentum sweep loses nothing it holds" begin
        # A field-free foil far from the momentum grid's ends: only the initial
        # Maxwellian's tails, below `cut` and flushed at the first sweep, cross
        # them (2e-21 of the electrons, 3e-18 of the ions). The spline's ringing
        # outside a column's support was dropped and counted as lost as well:
        # −1.9e-13 of the electrons and −8.6e-16 of the ions by 400 steps.
        vm = small_foil(; a0 = 0.0, scheme = SplineScheme())[1]
        N0 = [sum(density(vm, i))*vm.Δx for i in eachindex(vm.species)]
        for _ in 1:400
            step!(vm)
        end
        for (i, s) in enumerate(vm.species)
            @test 0 ≤ s.lost_p < 1e-16*N0[i]
            @test s.lost_p_energy ≥ 0
        end
    end

    @testset "the plasma and the probes are placed where the budget holds" begin
        vm, meta = small_foil(; a0 = 0.0)
        N, npml = meta.N, meta.npml
        build(; kw...) = VlasovMaxwell(; Δx = vm.Δx, Δt = vm.Δt, N, x_min = vm.x_min,
                                       pml = FDTD1D.PML(; N = npml, σ_max = 1.0, Δx = vm.Δx, Δt = vm.Δt),
                                       source = (y = (t, x) -> 0.0, z = (t, x) -> 0.0),
                                       species = Species[], kw...)
        # every face of the plasma an interior node, `npml + 2 : N − npml`: the
        # last one at `N − npml + 1`, the layer's first node, is refused, and
        # for the reason of the layer, not of the probes
        @test_throws "absorbing layer" build(; plasma_cells = (N - npml - 20):(N - npml))
        @test_throws "probes" build(; plasma_cells = (N - npml - 20):(N - npml - 1))
        @test_throws "absorbing layer" build(; plasma_cells = (npml + 1):(npml + 20))
        # the default left probe, `npml + 4`, inside a plasma that starts at the source
        @test_throws "probes" build(; plasma_cells = (npml + 2):(npml + 20))
        # probes past the source's cell, ahead of the plasma and behind it
        cells = (npml + 30):(N - npml - 30)
        @test build(; plasma_cells = cells) isa VlasovMaxwell
        @test_throws ArgumentError build(; plasma_cells = cells, probes = (npml + 3, N - npml - 2))
        @test_throws ArgumentError build(; plasma_cells = cells, probes = (npml + 40, N - npml - 2))
        @test_throws ArgumentError build(; plasma_cells = cells, probes = (npml + 4, N - npml - 40))
        @test_throws ArgumentError build(; plasma_cells = cells, probes = (npml + 4, N - npml + 1))
    end

    @testset "the step works in the data's type" begin
        function build(::Type{T}) where {T}
            Δx = T(2π/400); Δt = T(0.9)*Δx
            pml = FDTD1D.PML(; N = 50, σ_max = T(20), Δx, Δt)
            src = (y = (t, x) -> -T(0.5)*exp(-((t - 6)/3)^2)*sin(t - 6), z = (t, x) -> zero(T))
            e = momentum_edges(; pmin = -3.0, pmax = 3.0, Δp = 0.01, core = 0.3, ratio = 1.05, Δpmax = 0.1)
            fe = maxwellian(e, 0.002; n = 10.0)
            prof = [abs(c - 50) < 10 ? 1.0 : 0.0 for c in 1:100]
            sp = [Species(:e; charge = -1, mass = 1, edges = T.(e), f = T.(fe*prof'))]
            return VlasovMaxwell(; Δx, Δt, N = 600, x_min = zero(T), pml, source = src,
                                 plasma_cells = 250:349, species = sp, chunks = 1)
        end
        for T in (Float64, Float32)
            vm = build(T)
            N0 = sum(density(vm, 1))*vm.Δx
            peak = zero(T)
            for _ in 1:1500
                step!(vm)
                peak = max(peak, maximum(abs, vm.mesh.ey))
            end
            s = vm.species[1]
            drift = (sum(density(vm, 1))*vm.Δx + s.lost_x + s.lost_p + s.flushed - N0)/N0
            @test eltype(s.f) === T && eltype(vm.mesh.ey) === T && eltype(vm.ex) === T
            @test abs(drift) < 50eps(T)*length(s.pc)          # 2e-6 in Float32, 2e-15 in Float64
            @test peak > 0.3                                  # the pulse came through
        end
    end

    @testset "a Coulomb explosion conserves energy" begin
        # Electrons alone: a field of `n l` at the start, all of it turned into
        # kinetic energy as the foil explodes. No laser, so the budget closes on
        # the plasma's own two terms.
        # A window wide enough that no electron leaves it in the time run.
        vm, meta = laser_foil(; a0 = 0.0, ions = false, ξ = 0.5, chunks = 1, x_um = (-2.0, 2.0),
                              plasma_um = (-1.5, 1.6), pml_um = 0.4,
                              electron_grid = (pmin = -200.0, pmax = 200.0))
        E0 = kinetic_energy(vm) + longitudinal_energy(vm)
        for _ in 1:500
            step!(vm)
        end
        s = vm.species[1]
        E = kinetic_energy(vm) + longitudinal_energy(vm) + s.lost_x_energy + s.lost_p_energy
        println("  Coulomb explosion: field share at the start 1, now ",
                longitudinal_energy(vm)/E0, "; energy drift ", (E - E0)/E0)
        @test s.lost_x == 0 && s.lost_p == 0         # nothing has left the grid
        @test longitudinal_energy(vm) < 0.5*E0       # the explosion has happened
        @test abs(E - E0)/E0 < 1e-3
    end

    # The accuracy claims of `docs/src/verification.md`, "The relativistic laser
    # and foil": step 10's own setup, behind VASILEK_EXTENDED.
    @testset "extended: the laser and the foil" begin
        if get(ENV, "VASILEK_EXTENDED", "0") != "1"
            @info "laser-and-foil verification skipped; set VASILEK_EXTENDED=1 to run it"
            @test_skip "VASILEK_EXTENDED=1" == "set"
        else
            @testset "the pulse in vacuum" begin
                # Δx = 2 nm, 400 cells per wavelength, Courant number 0.937.
                r = vacuum_check()
                τp = r.meta.τp
                println("  vacuum: W_inc/formula − 1 = ", r.W_inc/r.W_formula - 1,
                        ", W(t₂)/W_inc − 1 = ", r.W_t2/r.W_inc - 1, ", G² leak ", r.G2_t2/r.W_inc,
                        ", peak ", r.amplitude, ", FWHM ", r.fwhm_fs, " fs, t_peak ", r.t_peak_fs,
                        " fs, layer ", r.R_layer)
                @test abs(r.W_inc/r.W_formula - 1) < 5e-4          # 1.3e-4
                @test abs(r.W_t2/r.W_inc - 1) < 1e-4               # 2.1e-6; README asks 5e-3
                @test r.G2_t2/r.W_inc < 1e-8                       # 8.5e-10
                # CEP 0 puts a node of the carrier on the envelope's peak: the
                # largest sample is a quarter cycle off it
                @test abs(r.amplitude - exp(-(π/2/τp)^2)) < 1e-4   # 0.99490 both
                @test abs(r.fwhm_fs - 11) < 0.01                   # 10.9991
                tf = 2τp*r.meta.u.T + 20/C_UM_PER_FS
                @test abs(r.t_peak_fs - tf) < 0.02                 # −0.007 fs, the source's 1.2Δx
                @test abs(r.phase_yee) < 1e-3                      # 2.0e-4 rad over 20 µm
                @test r.R_layer < 1e-10                            # 2.3e-13
            end

            @testset "linear reflection and transmission" begin
                # a₀ = 1e-3 on 10, 30 and 90 nm at 40 n_c, ions free, 1 keV. The
                # foil expands for 85 fs before the pulse's peak, the thinnest by
                # more than its thickness, so the reference is Fresnel's layer
                # stack of the run's own electron profile, averaged over the pulse.
                # Measured: R +0.16%, +0.24%, +0.25%; T +0.36%, +0.49%, +1.0%.
                for l in (10.0, 30.0, 90.0)
                    q = fresnel_check(l)
                    println("  Fresnel, $l nm: R = ", q.R, " against ", q.profile[1], ", T = ", q.T,
                            " against ", q.profile[2], ", R + T − 1 = ", q.R + q.T - 1,
                            ", front at ", q.v_front, " c")
                    @test abs(q.R/q.profile[1] - 1) < 4e-3
                    @test abs(q.T/q.profile[2] - 1) < 2e-2
                    # a surface coming at the pulse returns more than it got:
                    # R + T − 1 is positive and under the Doppler bound 2|v|
                    @test 0 ≤ q.R + q.T - 1 ≤ 2*abs(q.v_front)
                    @test abs(q.areal - 1) < 2e-3
                end
            end

            @testset "the nonlinear budget closes" begin
                # a₀ = 20, ξ = 0.972 at 4 nm: balance 0.99693, residual −3.1e-3
                # (2 nm: 0.99870, −1.3e-3).
                n = foil_transmission(; dx_nm = 4.0)
                println("  a₀ = 20, ξ = 0.972, 4 nm: T = ", n.T, " R = ", n.R, ", balance ", n.balance,
                        ", residual ", n.residual)
                @test abs(n.balance - 1) < 5e-3
                @test abs(n.residual) < 5e-3
                @test n.vm.species[1].lost_p == 0
                @test 0.7 < n.T < 0.8                       # 0.7487
            end
        end
    end
end
