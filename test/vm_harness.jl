# The laser-and-foil problem for `VlasovMaxwell1D1P`, shared by the tests, the
# verification studies in `verification/` and anything else that runs it, so
# that a number in `docs/src/verification.md` and the test that pins it come
# from one setup rather than two copies.
#
# Units: the module's (time 1/ω₀, length c/ω₀, fields mₑcω₀/e, densities n_c);
# the arguments here are in the laboratory's (µm, fs, nm, eV), converted once.

using Vasilek
using Vasilek.VlasovMaxwell1D1P
using Vasilek: FDTD1D

const VM_HARNESS = true

"""
    LaserUnits(λ_µm = 0.8)

The conversions between the module's units and the laboratory's for a laser
of wavelength `λ_µm`: `L` µm per `c/ω₀`, `T` fs per `1/ω₀`.
"""
struct LaserUnits
    λ::Float64
    L::Float64
    T::Float64
end
const C_UM_PER_FS = 0.299792458
LaserUnits(λ = 0.8) = LaserUnits(λ, λ/2π, λ/2π/C_UM_PER_FS)

"Electron rest energy, eV."
const MEC2_EV = 510998.95
"`12 amu − 6 mₑ` in electron masses, the C⁶⁺ mass of `spawn_cs.py`."
const C6_MASS = 12*1822.888486209 - 6

"""
    gaussian_pulse(t; a0, τp, cep = 0)

The laser of step 10 (`README.md` §2), in the module's units: the field that
crosses the entry plane at time `t`,

    a0·exp(−((t − 2τp)/τp)²)·sin(t − 2τp + cep),     zero before t = 0.
"""
gaussian_pulse(t; a0, τp, cep = 0.0) =
    t < 0 ? zero(t) : a0*exp(-((t - 2τp)/τp)^2)*sin(t - 2τp + cep)

"""
    laser_foil(; a0 = 20, ξ = 0.972, l_star_nm = 90, n0 = 40, Te_eV = 1000,
               Ti_eV = Te_eV, dx_nm = 2, steps_per_tenth_fs = nothing, cfl = 0.937,
               x_um = (-20, 15), plasma_um = (-5, 10), pml_um = 1.0, pml_R = 1e-12,
               τ_fs = 11, λ_µm = 0.8, cep = 0, ions = true, Zi = 6, mi = C6_MASS,
               electron_grid = (;), ion_grid = (;), scheme = PFCScheme(),
               chunks = Threads.nthreads(), cut = 1e-13, vacuum = false,
               l_nm = ξ*l_star_nm)

A run of step 10: a plane laser pulse entering at `x_um[1]`, a foil of
thickness `l_nm` at `x ∈ [0, l]` with electrons at `n0` and fully ionised ions
of charge `Zi` and mass `mi` at `n0/Zi`, both Maxwellian, `f` kept on
`plasma_um`, absorbing layers `pml_um` thick outside `x_um`.

The foil's areal density is exact: a cell the surface cuts holds the fraction
of it inside. `Δt` divides 0.1 fs exactly (`steps_per_tenth_fs`, by default the
fewest steps with `cfl ≤ 0.937`) so that every output time of the study is a
step. Returns `(vm, meta)`; `meta` carries the units, the indices of the
nodes and cells at given positions and what was set.
"""
function laser_foil(; a0 = 20.0, ξ = 0.972, l_star_nm = 90.0, n0 = 40.0,
                    Te_eV = 1000.0, Ti_eV = Te_eV, dx_nm = 2.0, steps_per_tenth_fs = nothing,
                    cfl = 0.937, x_um = (-20.0, 15.0), plasma_um = (-5.0, 10.0),
                    pml_um = 1.0, pml_R = 1e-12, τ_fs = 11.0, λ_µm = 0.8, cep = 0.0,
                    ions = true, Zi = 6, mi = C6_MASS, electron_grid = (;), ion_grid = (;),
                    scheme = PFCScheme(), chunks = 4Threads.nthreads(), cut = 1e-13,
                    vacuum = false, l_nm = ξ*l_star_nm)
    u = LaserUnits(λ_µm)
    Δx = dx_nm*1e-3/u.L
    # Δt: a whole number of steps per 0.1 fs, Courant number at most `cfl`
    tenth = 0.1/u.T
    nsub = steps_per_tenth_fs === nothing ? ceil(Int, tenth/(cfl*Δx)) : steps_per_tenth_fs
    Δt = tenth/nsub
    # the grid: node `ks = npml + 2` at x_um[1], the source's node
    npml = round(Int, pml_um*1e3/dx_nm)
    ks = npml + 2
    node_at(x) = ks + round(Int, (x - x_um[1])*1e3/dx_nm)
    kend = node_at(x_um[2])
    N = kend + npml - 1          # the right layer's nodes kend+1 … N, its wall N+1
    x_min = x_um[1]/u.L - (ks - 1)*Δx
    k1, k2 = node_at(plasma_um[1]), node_at(plasma_um[2])
    plasma_cells = k1:(k2 - 1)
    Nc = length(plasma_cells)
    # the layer: σ_max for a reflection of `pml_R` at normal incidence through
    # a cubic profile, `R = exp(−2∫σ dx) = exp(−σ_max L/2)`
    L = npml*Δx
    σ_max = 2*log(1/pml_R)/L
    pml = FDTD1D.PML(; N = npml, σ_max, Δx, Δt)
    τp = τ_fs/sqrt(2log(2))/u.T
    # `Yee1D` launches minus its source
    source = (y = (t, x) -> -gaussian_pulse(t; a0, τp, cep), z = (t, x) -> 0.0)
    # the foil's fraction in each plasma cell
    xc_face = [x_um[1]/u.L + (k - ks)*Δx for k in k1:k2]     # faces, module units
    l = l_nm*1e-3/u.L
    frac = [max(0.0, min(xc_face[c+1], l) - max(xc_face[c], 0.0))/Δx for c in 1:Nc]
    species = Species[]
    if !vacuum
        eg = merge((pmin = -80.0, pmax = 300.0, Δp = 0.01, core = 0.5, ratio = 1.05,
                    Δpmax = 1.0), electron_grid)
        ee = momentum_edges(; eg...)
        fe = maxwellian(ee, Te_eV/MEC2_EV; n = n0)
        push!(species, Species(:electron; charge = -1, mass = 1, edges = ee, f = fe*frac'))
        if ions
            θi = Ti_eV/(MEC2_EV*mi)
            pth = sqrt(θi)
            ig = merge((pmin = -0.03, pmax = 0.08, Δp = pth/4, core = 10pth, ratio = 1.05,
                        Δpmax = 1e-3), ion_grid)
            ei = momentum_edges(; ig...)
            fi = maxwellian(ei, θi; n = n0/Zi)
            push!(species, Species(:ion; charge = Zi, mass = mi, edges = ei, f = fi*frac'))
        end
    end
    vm = VlasovMaxwell(; Δx, Δt, N, x_min, pml, source, plasma_cells, species, scheme,
                       chunks, cut, probes = (ks + 2, kend - 2))
    xs = x_min .+ Δx .* (0:N)
    meta = (; u, Δx, Δt, N, ks, kend, npml, σ_max, k1, k2, plasma_cells, τp, a0, ξ, l_nm,
            n0, Te_eV, Ti_eV, dx_nm, nsub, x_um, plasma_um, Zi, mi, xs, node_at,
            S = n0/a0, cfl = Δt/Δx)
    return vm, meta
end

"""
    incident_energy(a0, τp)

`∫F²dt` of [`gaussian_pulse`](@ref) from `t = 0`, in the module's units: the
energy per area it carries across the entry plane, `n_c mₑc² (c/ω₀)`.
"""
function incident_energy(a0, τp; cep = 0.0)
    t = range(0, 8τp; length = 400_001)
    y = gaussian_pulse.(t; a0, τp, cep) .^ 2
    return (sum(y) - (y[1] + y[end])/2)*step(t)
end

"""
    waves(vm)

`F = (E_y + B_z)/2` and `G = (E_y − B_z)/2` at every node at the step's end,
`B_z` the mean of the node's two cells at the two half steps.
"""
function waves(vm)
    ey, hz, hp = vm.mesh.ey, vm.mesh.hz, vm.hprev
    N = vm.mesh.N
    B = zeros(N + 1)
    for k in 2:N
        B[k] = (hz[k-1] + hz[k] + hp[k-1] + hp[k])/4
    end
    return (ey .+ B)./2, (ey .- B)./2
end

"""
    window_integral(y, meta, x1_um, x2_um)

`Σ y Δx` over the nodes in `[x1_um, x2_um]`, in the module's units.
"""
function window_integral(y, meta, x1, x2)
    k1, k2 = meta.node_at(x1), meta.node_at(x2)
    return sum(view(y, k1:k2))*meta.Δx
end

"""
    steps_to(meta, t_fs)

The number of steps from `t = 0` to `t_fs`, which must be a step.
"""
function steps_to(meta, t_fs)
    n = t_fs/0.1*meta.nsub
    isapprox(n, round(n); atol = 1e-6) || throw(ArgumentError("$t_fs fs is not a step"))
    return round(Int, n)
end

"""
    vacuum_check(; a0 = 20, dx_nm = 2, t1 = 60, t2 = 110.4, t_late = 175, kw...)

The laser of [`laser_foil`](@ref) without a foil. Returns, in the module's
units unless named otherwise:

  * `W_formula`: `∫F²dt` of the formula at the entry plane, [`incident_energy`](@ref);
  * `W_inc`: `∫F²dx` over `[−20, −1]` µm at `t1` (README §5.1), and `W_t2` over
    `[1, 15]` µm at `t2`, the same pulse after 50.4 fs of flight;
  * `G2_t2`: `∫G²dx` over the whole interior at `t2`, the source's leak to the left;
  * `amplitude`: the peak of `|E_y|` at the node `x = 0` against `a0`, and
    `fwhm_fs`, the intensity envelope's full width there (the formula: 11 fs);
  * `t_peak_fs`: when the envelope at `x = 0` peaks (the formula: `t_f = 85.4` fs);
  * `phase`: the carrier's phase at `x = 0` against the formula's delayed by
    `20 µm/c`, fitted, rad, and `phase_yee`, what Yee's dispersion predicts;
  * `R_layer`: `∫G²dx` over the interior at `t_late`, once the pulse has
    entered the right absorbing layer, over `W_inc`.
"""
function vacuum_check(; a0 = 20.0, dx_nm = 2.0, t1 = 60.0, t2 = 110.4, t_late = 175.0, kw...)
    vm, meta = laser_foil(; a0, dx_nm, vacuum = true, kw...)
    u = meta.u
    k0 = meta.node_at(0.0)
    ts = Float64[]; e0 = Float64[]
    W_inc = W_t2 = G2_t2 = NaN
    n1, n2, n3 = steps_to(meta, t1), steps_to(meta, t2), steps_to(meta, t_late)
    for n in 1:n3
        step!(vm)
        push!(ts, vm.t); push!(e0, vm.mesh.ey[k0])
        if n == n1 || n == n2
            F, G = waves(vm)
            if n == n1
                W_inc = window_integral(F.^2, meta, -20.0, -1.0)
            else
                W_t2 = window_integral(F.^2, meta, 1.0, 15.0)
                G2_t2 = window_integral(G.^2, meta, meta.x_um...)
            end
        end
    end
    F, G = waves(vm)
    R_layer = window_integral(G.^2, meta, meta.x_um...)/W_inc
    # the envelope at x = 0 through the analytic signal of the carrier
    env = _envelope(e0)
    i = argmax(env)
    half = env[i]^2/2
    lo = findlast(<(half), view(env .^ 2, 1:i))
    hi = i - 1 + findfirst(<(half), view(env .^ 2, i:length(env)))
    cross(j) = ts[j] + (half - env[j]^2)/(env[j+1]^2 - env[j]^2)*(ts[j+1] - ts[j])
    fwhm = (cross(hi - 1) - cross(lo))*u.T
    t_peak = _parabola_peak(ts, env, i)*u.T
    # the carrier's phase against the formula delayed by the flight
    L = 20.0/u.L
    ref = [gaussian_pulse(t - L; a0, τp = meta.τp) for t in ts]
    quad = [gaussian_pulse(t - L; a0, τp = meta.τp, cep = π/2) for t in ts]
    phase = atan(sum(e0 .* quad), sum(e0 .* ref))
    # Yee's phase delay over L: the grid's wavenumber at ω = 1 against 1
    Δx, Δt = meta.Δx, meta.Δt
    kyee = 2/Δx*asin(Δx/Δt*sin(Δt/2))
    phase_yee = -(kyee - 1)*L
    return (; W_formula = incident_energy(a0, meta.τp), W_inc, W_t2, G2_t2,
            amplitude = maximum(abs, e0)/a0, fwhm_fs = fwhm, t_peak_fs = t_peak,
            phase, phase_yee, R_layer, meta)
end

# The envelope of a sampled carrier at ω = 1, |analytic signal|, through the FFT
# of the real signal with the negative frequencies removed.
function _envelope(y)
    n = length(y)
    Y = Vasilek.PoissonFourier1D.FFTW.fft(y)
    h = zeros(n); h[1] = 1
    if iseven(n)
        h[n ÷ 2 + 1] = 1; h[2:n÷2] .= 2
    else
        h[2:(n+1)÷2] .= 2
    end
    return abs.(Vasilek.PoissonFourier1D.FFTW.ifft(Y .* h))
end

# The vertex of the parabola through the three samples around `i`.
function _parabola_peak(t, y, i)
    y₋, y₀, y₊ = y[i-1], y[i], y[i+1]
    d = (y₋ - y₊)/(2*(y₋ - 2y₀ + y₊))
    return t[i] + d*(t[i+1] - t[i])
end

"""
    fresnel_slab(ω, n, l)

Reflectance and transmittance `(R, T)` of a uniform cold slab of density `n`
(in `n_c` at `ω₀ = 1`) and thickness `l` (`c/ω₀`) in vacuum, at normal
incidence and frequency `ω`: `ε = 1 − n/ω²`, through the characteristic matrix
of the layer, `δ = √ε ωl`.
"""
function fresnel_slab(ω, n, l)
    N = sqrt(complex(1 - n/ω^2))
    δ = N*ω*l
    m11 = m22 = cos(δ)
    m12 = -im*sin(δ)/N
    m21 = -im*N*sin(δ)
    d = m11 + m12 + m21 + m22
    r = (m11 + m12 - m21 - m22)/d
    t = 2/d
    return abs2(r), abs2(t)
end

"""
    fresnel_pulse(n, l, τp)

[`fresnel_slab`](@ref) averaged over the power spectrum of
[`gaussian_pulse`](@ref), `∝ exp(−(ω − 1)²τp²/2)`: what the energies of the
reflected and transmitted pulses give.
"""
function fresnel_pulse(n, l, τp)
    ωs = range(1 - 8/τp, 1 + 8/τp; length = 4001)
    w = @. exp(-(ωs - 1)^2*τp^2/2)
    RT = fresnel_slab.(ωs, n, l)
    return sum(w .* first.(RT))/sum(w), sum(w .* last.(RT))/sum(w)
end

"""
    thin_foil_R(n, l)

The thin-foil limit `R = ζ²/(1 + ζ²)`, `ζ = n l/2` (`π n l/(n_c λ)` in the
laboratory's units).
"""
thin_foil_R(n, l) = (ζ = n*l/2; ζ^2/(1 + ζ^2))

"""
    fresnel_profile(ω, n, h)

`(R, T)` of a stack of cold layers of densities `n` (in `n_c` at `ω₀ = 1`),
each `h` thick, at frequency `ω`: the product of their characteristic
matrices. A uniform slab of `m` layers is [`fresnel_slab`](@ref).
"""
function fresnel_profile(ω, n, h)
    m11, m12, m21, m22 = complex(1.0), complex(0.0), complex(0.0), complex(1.0)
    for nk in n
        N = sqrt(complex(1 - nk/ω^2))
        δ = N*ω*h
        a11, a12, a21, a22 = cos(δ), -im*sin(δ)/N, -im*N*sin(δ), cos(δ)
        m11, m12, m21, m22 = m11*a11 + m12*a21, m11*a12 + m12*a22,
                             m21*a11 + m22*a21, m21*a12 + m22*a22
    end
    d = m11 + m12 + m21 + m22
    return abs2((m11 + m12 - m21 - m22)/d), abs2(2/d)
end

"""
    fresnel_check(l_nm; a0 = 1e-3, dx_nm = 2, t2 = 110.4, plasma_um = (-0.5, 0.6), kw...)

The linear reflection and transmission of a foil `l_nm` thick at `n0 = 40`,
measured as step 10 measures them, `R = ∫G²dx` over `[−20, −1]` µm and
`T = ∫F²dx` over `[1, 15]` µm at `t2` over the vacuum run's `W_inc`, against

  * `slab`: [`fresnel_slab`](@ref) of the foil as it was set up, averaged over
    the pulse's spectrum ([`fresnel_pulse`](@ref));
  * `thin`: the thin-foil limit [`thin_foil_R`](@ref);
  * `profile`: the layer stack of the electrons' own profile
    ([`fresnel_profile`](@ref)), averaged over the spectrum and over the
    profiles every 2 fs around the pulse's arrival `t_f`, weighted by its
    intensity there, `exp(−2((t − t_f)/τ_p)²)`: what the foil, expanding as the
    pulse crosses it, should do if it stood still each moment;
  * `v_front`: the speed, in `c`, of the point where `nₑ` first reaches `n_c`
    from the left, fitted over `t_f ± τ_p`: a surface coming at the pulse
    returns it Doppler-shifted, `R·(1 + 2v)` to first order.

The electrons' momentum grid is cut to `±2`, which a quiver of `1e-3` and a
temperature of 1 keV stay well inside.
"""
function fresnel_check(l_nm; a0 = 1e-3, dx_nm = 2.0, t2 = 110.4, plasma_um = (-0.5, 0.6), kw...)
    vac, mv = laser_foil(; a0, dx_nm, vacuum = true,
                         (k => v for (k, v) in pairs(kw) if k in (:steps_per_tenth_fs, :cfl))...)
    for _ in 1:steps_to(mv, 60.0); step!(vac); end
    F, _ = waves(vac)
    W_inc = window_integral(F.^2, mv, -20.0, -1.0)
    vm, meta = laser_foil(; a0, dx_nm, l_nm, ξ = l_nm/90, plasma_um,
                          electron_grid = (pmin = -2.0, pmax = 2.0), kw...)
    u = meta.u
    τp_fs = meta.τp*u.T
    tf = round(2τp_fs + 20/C_UM_PER_FS; digits = 1)
    times = [round(tf + 2k; digits = 1) for k in -8:8]
    marks = Dict(steps_to(meta, t) => t for t in times)
    profiles = Dict{Float64, Vector{Float64}}()
    for n in 1:steps_to(meta, t2)
        step!(vm)
        haskey(marks, n) && (profiles[marks[n]] = density(vm, 1))
    end
    F, G = waves(vm)
    T = window_integral(F.^2, meta, 1.0, 15.0)/W_inc
    R = window_integral(G.^2, meta, -20.0, -1.0)/W_inc
    l = l_nm*1e-3/u.L
    slab = fresnel_pulse(meta.n0, l, meta.τp)
    ωs = range(1 - 8/meta.τp, 1 + 8/meta.τp; length = 401)
    w = @. exp(-(ωs - 1)^2*meta.τp^2/2)
    Rp = Tp = W = 0.0
    xc = collect(cell_x(vm)) .* u.L
    fronts = Tuple{Float64, Float64}[]
    for t in times
        ne = profiles[t]
        RT = [fresnel_profile(ω, ne, meta.Δx) for ω in ωs]
        wt = exp(-2((t - tf)/τp_fs)^2)
        Rp += wt*sum(w .* first.(RT))/sum(w)
        Tp += wt*sum(w .* last.(RT))/sum(w)
        W += wt
        if abs(t - tf) ≤ τp_fs
            i = findfirst(≥(1.0), ne)
            # the crossing between cells i-1 and i, linearly
            xf = xc[i-1] + (1 - ne[i-1])/(ne[i] - ne[i-1])*(xc[i] - xc[i-1])
            push!(fronts, (t, xf))
        end
    end
    ts, xs = first.(fronts), last.(fronts)
    v_front = sum((ts .- sum(ts)/length(ts)) .* (xs .- sum(xs)/length(xs)))/
              sum((ts .- sum(ts)/length(ts)).^2)/C_UM_PER_FS
    ne_f = profiles[tf]
    return (; R, T, A = 1 - R - T, slab, thin = thin_foil_R(meta.n0, l), profile = (Rp/W, Tp/W),
            v_front, W_inc, meta, areal = sum(ne_f)*meta.Δx/(meta.n0*l))
end

"""
    foil_transmission(; a0 = 20, ξ = 0.972, dx_nm = 2, t2 = 110.4, W_inc = nothing, kw...)

The nonlinear run of step 10 to `t2`: `T`, `R` and `A = 1 − T − R` measured as
[`fresnel_check`](@ref) measures them, and the energy budget of
[`energy_account`](@ref) over `W_inc`: `balance`, the field, kinetic, lost and
outgoing energy together, and `residual`, their sum less what came in through
the left probe. `W_inc` is taken from a vacuum run on the same grid unless
given. Keywords go to [`laser_foil`](@ref).
"""
function foil_transmission(; a0 = 20.0, ξ = 0.972, dx_nm = 2.0, t2 = 110.4, W_inc = nothing,
                           kw...)
    if W_inc === nothing
        vac, mv = laser_foil(; a0, dx_nm, vacuum = true,
                             (k => v for (k, v) in pairs(kw) if k in (:steps_per_tenth_fs, :cfl))...)
        for _ in 1:steps_to(mv, 60.0); step!(vac); end
        F, _ = waves(vac)
        W_inc = window_integral(F.^2, mv, -20.0, -1.0)
    end
    vm, meta = laser_foil(; a0, ξ, dx_nm, kw...)
    for _ in 1:steps_to(meta, t2)
        step!(vm)
    end
    F, G = waves(vm)
    T = window_integral(F.^2, meta, 1.0, 15.0)/W_inc
    R = window_integral(G.^2, meta, -20.0, -1.0)/W_inc
    acc = energy_account(vm)
    balance = (acc.field + acc.longitudinal + acc.kinetic + acc.lost + acc.out_left +
               acc.out_right)/W_inc
    return (; T, R, A = 1 - T - R, balance, residual = acc.residual/W_inc, acc, W_inc, vm, meta)
end
