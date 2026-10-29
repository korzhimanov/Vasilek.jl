"""
    VlasovMaxwell1D1P

Relativistic Vlasov–Maxwell in one dimension of space and one of momentum, for
any number of species, with a linearly polarised laser through
[`FDTD1D`](@ref Vasilek.FDTD1D).

The model is that of `korzhimanov/vlasov-solver`, the C++ code this module
ports. Nothing depends on `y` or `z`, so each particle keeps its transverse
canonical momentum, zero ahead of the pulse: `p⊥ + (q/e)(mₑ/m)a = 0`, with `a`
the normalised vector potential `eA_y/(mₑc)`. A species of charge `q` (in `e`)
and mass `m` (in `mₑ`) then has, in its own momentum unit `mc`,

    γ = √(1 + p² + κ²a²),     κ = q/m,
    ∂f/∂t + ∂(v f)/∂x + ∂(F f)/∂p = 0,     v = ∂γ/∂p,  F = κEₓ − ∂γ/∂x,

and the fields obey Yee's update for `(E_y, B_z)` with the transverse current
`j_y = −Σ q κ a ∫ f/γ dp`, `∂a/∂t = −E_y`, and Ampère's law `∂Eₓ/∂t = −jₓ`.
Units are those of the electromagnetic section of `docs/src/normalization.md`:
time `1/ω₀`, length `c/ω₀`, fields `mₑcω₀/e`, densities `n_c`.

# The discretisation

* `f[j, c]` is the average over momentum cell `j` and space cell `c`. The space
  cells are the Yee cells: their faces are the electric nodes, where `a`,
  `E_y` and `Eₓ` live.
* The Hamiltonian `γ − 1` is evaluated at the cell corners, the momentum
  **edges** of each face, and `v` and `F` are its differences: `v` across a
  momentum cell along a face, `F` across a space cell along an edge. The
  velocity field is then divergence-free on the grid, exactly -- the two
  differences cancel cell by cell -- as the Hamiltonian flow it stands for is.
* Each direction is a flux-form, conservative sweep: in `x` the PFC flux of
  Filbet, Sonnendrücker and Bertrand (2001), the displacement being under half a
  cell; in `p`, where one step moves `f` by tens of the finest cells, the same
  reconstruction integrated from the foot of each edge's characteristic, a
  conservative semi-Lagrangian step with no Courant limit
  ([`psweep!`](@ref)). Both are positive, and conserve the number of particles
  to round-off apart from what crosses the ends of the grid, which is counted.
* The longitudinal field is advanced by the charge each `x` sweep moves across
  each face, so that Gauss's law holds at every step to round-off if it held at
  the start, charge leaving the plasma's grid included: it stays where it left.
* A step is `X(Δt/2) P(Δt) X(Δt/2)` with `a` frozen at the half step,
  `a^{n+1/2} = a^{n-1/2} − Δt E_y^n`, which is second order in time. The
  transverse current is the mean of the currents before and after `P`.

# Splitting a relativistic flow

Huot, Ghizzo, Bertrand, Sonnendrücker and Coulaud (J. Comput. Phys. 185, 512,
2003) found the time splitting of this system unstable: with `γ` depending on
both `x` and `p`, neither split flow is incompressible, and a split step that
carries `f` along each sub-flow's characteristics as a value, as a backward
semi-Lagrangian interpolation does, conserves neither the particles nor the
phase-space volume. The C++ code splits all the same, `P(Δt)` then `X(Δt)`,
but in conservative form, each cell's content moved and deposited, so that the
particles are conserved by each sub-step whatever the compression. This module
keeps that, makes the two sub-flows' compressions cancel exactly on the grid
(the corner Hamiltonian above), and symmetrises the splitting; the nonlinear
runs of the transparency study close their energy budget to 0.13% at 2 nm.

# Supports

Every column of `f` carries the momentum cells outside which it is zero, and
the sweeps work on those alone: away from the foil a column holds a beam a
few cells wide out of hundreds. So that the supports mean something, values
below `cut` times a species' initial maximum are set to zero after each
momentum sweep, and their mass counted; with `cut = 1e-13` it was 3e-7 of the
electrons in the transparency study's reference run.

`step!` allocates nothing on one thread. With more, it spreads each sweep over
the threads in contiguous chunks, `4 Threads.nthreads()` by default and taken
dynamically; the result does not depend on their number.
"""
module VlasovMaxwell1D1P

using ..FDTD1D: YeeMesh1D, Yee1D, PML, advance!
using ..Advection: slope_limit

export Species, momentum_edges, maxwellian, VlasovMaxwell, PFCScheme, SplineScheme,
       step!, density, density!, distribution, support, kinetic_energy, field_energy,
       longitudinal_energy, node_x, cell_x, plasma_nodes, energy_account

# ------------------------------------------------------------------ grids

"""
    momentum_edges(; pmin, pmax, Δp, core = 0, ratio = 1.05, Δpmax = Inf)

Edges of a momentum grid: cells of width `Δp` over `[-core, core]` (rounded to
a whole number of cells, at least one either side of an edge at 0), then cells
growing outward by `ratio` each, never wider than `Δpmax`, to the first edge at
or past `pmin` and `pmax`.

Growth by `ratio` keeps the relative resolution `Δp/|p|` near `ratio − 1` far
out, which is what an energy spectrum needs; `Δpmax` caps it where the
dynamics at large `|p|` still wants absolute resolution.
"""
function momentum_edges(; pmin, pmax, Δp, core = 0.0, ratio = 1.05, Δpmax = Inf)
    pmin < 0 < pmax || throw(ArgumentError("need pmin < 0 < pmax, got [$pmin, $pmax]"))
    Δp > 0 && ratio ≥ 1 && Δpmax ≥ Δp || throw(ArgumentError(
        "need Δp > 0, ratio ≥ 1 and Δpmax ≥ Δp; got $Δp, $ratio, $Δpmax"))
    n = max(1, round(Int, core/Δp))
    side(bound) = begin
        e = [k*Δp for k in 0:n]
        w = float(Δp)
        while e[end] < bound
            w = min(w*ratio, Δpmax)
            push!(e, e[end] + w)
        end
        e
    end
    right, left = side(pmax), side(-pmin)
    return vcat(-reverse(left[2:end]), right)
end

"""
    maxwellian(edges, θ; n = 1)

Cell averages over `edges` of the relativistic one-dimensional equilibrium
`∝ exp(−(√(1 + p²) − 1)/θ)`, `θ = T/(mc²)`, scaled so that `Σ f Δp = n` on the
grid exactly. It is a function of `γ` alone, so it is a steady state of the
field-free system, and for `θ ≪ 1` the Gaussian of temperature `θ`.
"""
function maxwellian(edges, θ; n = 1.0)
    θ > 0 || throw(ArgumentError("θ = $θ: the temperature must be positive"))
    # Gauss–Legendre with five nodes on each cell, enough for a Gaussian
    # resolved by a few cells and exact to round-off for one resolved by more.
    xg = (-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
          0.9061798459386640)
    wg = (0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
          0.4786286704993665, 0.2369268850561891)
    np = length(edges) - 1
    f = zeros(np)
    for j in 1:np
        a, b = edges[j], edges[j+1]
        m, h = (a + b)/2, (b - a)/2
        s = 0.0
        for (x, w) in zip(xg, wg)
            p = m + h*x
            s += w*exp(-(p*p/(1 + sqrt(1 + p*p)))/θ)
        end
        f[j] = s/2           # the average over the cell
    end
    Δp = diff(edges)
    f .*= n/sum(f .* Δp)
    return f
end

"""
    Species(name; charge, mass, edges, f)

A species for [`VlasovMaxwell`](@ref): its charge in `e` (−1 for electrons),
mass in `mₑ`, momentum grid `edges` in its own `mc`, and the initial
distribution `f`, a matrix `[momentum cell, space cell]` of cell averages in
`n_c` per unit momentum, over the plasma's space cells. `f ≥ 0`.
"""
struct Species{T<:AbstractFloat}
    name::Symbol
    charge::T
    mass::T
    edges::Vector{T}
    f::Matrix{T}
    function Species(name; charge, mass, edges, f)
        T = float(promote_type(typeof(charge), typeof(mass), eltype(edges), eltype(f)))
        mass > 0 || throw(ArgumentError("species $name: mass must be positive, got $mass"))
        length(edges) ≥ 4 && issorted(edges; lt = ≤) || throw(ArgumentError(
            "species $name: edges must be strictly increasing, at least three cells"))
        size(f, 1) == length(edges) - 1 || throw(DimensionMismatch(
            "species $name: f has $(size(f, 1)) momentum cells, edges $(length(edges) - 1)"))
        all(≥(0), f) || throw(DomainError(minimum(f), "species $name: f must be ≥ 0"))
        return new{T}(Symbol(name), T(charge), T(mass), Vector{T}(edges), Matrix{T}(f))
    end
end

# ------------------------------------------------------------------ schemes

"""
    PFCScheme()

The positive, flux-conservative reconstruction in both directions (default).
"""
struct PFCScheme end

"""
    SplineScheme()

The cubic-spline semi-Lagrangian step in both directions, in conservative
form: the primitive `∫f` is interpolated by a cubic spline with zero slope at
the ends and differenced between the feet of the characteristics. Third order,
not positive. For comparison with [`PFCScheme`](@ref).
"""
struct SplineScheme end

# ------------------------------------------------------------------ state

"""
    SpeciesState

A species as [`VlasovMaxwell`](@ref) steps it. `f` and its scratch `g` are
`[momentum cell, space column]` with two empty ghost columns at either end;
space cell `c` is column `c + 2`.

Every column carries its **support**, `jlo[col]:jhi[col]`, the momentum cells
outside which it is zero (`jlo > jhi` for an empty one), and `g` its own. A
plasma column away from the foil holds a beam a few cells wide out of
hundreds, and the sweeps work on the supports alone. `lo:hi` are the cells
whose support may be non-empty, `glo:ghi` those of `g`.

`A2` is `(κa)²` on every face for the current step, from which `γ − 1` at any
momentum is one square root; `w` holds each cell's two sums `Σ f Δp/γ`
towards its faces, the current's raw material, for the cells `wlo:whi`. The
accumulators count, in `n_c·(c/ω₀)` and `n_c mₑc²·(c/ω₀)`, what left through
the ends of `x` and of `p`, and what was flushed below `cut`.
"""
mutable struct SpeciesState{T<:AbstractFloat}
    const name::Symbol
    const q::T
    const m::T
    const κ::T
    const pe::Vector{T}
    const pc::Vector{T}
    const Δp::Vector{T}
    const iΔp::Vector{T}
    const pe2::Vector{T}
    const pc2::Vector{T}
    const ξ::Vector{T}
    const r3::Vector{T}      # 1/(d₋ + d₀ + d₊) per momentum cell
    const r2p::Vector{T}     # 1/(d₀ + d₊)
    const r2m::Vector{T}     # 1/(d₋ + d₀)
    const j0::Int            # the first momentum cell above the edge at 0
    const fmax::T
    const cut::T
    f::Matrix{T}
    g::Matrix{T}
    jlo::Vector{Int}
    jhi::Vector{Int}
    gjlo::Vector{Int}
    gjhi::Vector{Int}
    const A2::Vector{T}
    const w::Matrix{T}
    lo::Int
    hi::Int
    glo::Int
    ghi::Int
    wlo::Int
    whi::Int
    lost_x::T
    lost_x_energy::T
    lost_p::T
    lost_p_energy::T
    flushed::T
end

struct ChunkBuffers{T<:AbstractFloat}
    Φa::Vector{T}
    Φb::Vector{T}
    C::Vector{T}
    M::Vector{T}
    y2::Vector{T}
    tmp::Vector{T}
    tmp2::Vector{T}
    xpiv::Vector{T}
    jx::Vector{T}
    acc::Vector{T}
    rng::Vector{Int}     # the rows Φa and Φb hold: (alo, ahi, blo, bhi)
end

function ChunkBuffers{T}(np, nf) where {T}
    L = max(np, nf) + 5
    return ChunkBuffers{T}(zeros(T, np + 1), zeros(T, np + 1), zeros(T, L), zeros(T, L),
                           zeros(T, L), zeros(T, L), zeros(T, L), zeros(T, L), zeros(T, nf),
                           zeros(T, 6), [1, 0, 1, 0])
end

"""
    VlasovMaxwell(; Δx, Δt, N, x_min, pml, source, plasma_cells, species,
                  scheme = PFCScheme(), chunks = 4Threads.nthreads(), cut = 1e-13,
                  probes = (pml.N + 4, N - pml.N - 2))

The state of a run, at `t = 0`, ready for [`step!`](@ref).

  * `Δx`, `Δt`, `N`, `x_min`, `pml`, `source`: the Yee grid of `N` cells, its
    first node at `x_min`, and its operator, [`Yee1D`](@ref
    Vasilek.FDTD1D.Yee1D); the source is injected one-way at node `pml.N + 2`,
    and launches a wave of **minus** its value.
  * `plasma_cells`: the Yee cells `f` lives on, inside both absorbing layers.
    What crosses their outer faces is removed and counted.
  * `species`: [`Species`](@ref) over those cells.
  * `cut`: values of `f` below `cut` times the species' initial maximum are
    set to zero after each momentum sweep, and counted, so that the cells a
    run works on are those the plasma has reached rather than those its
    exponentially small tails have.
  * `probes`: the two nodes the energy flux is integrated at, for
    [`energy_account`](@ref).

`E_y`, `B_z`, `a` and `Eₓ` start at zero, `Eₓ` then integrated from the
plasma's charge, so a non-neutral start has its field.
"""
mutable struct VlasovMaxwell{T<:AbstractFloat, OP<:Yee1D, S}
    const mesh::YeeMesh1D{T}
    const op::OP
    const scheme::S
    const Δx::T
    const Δt::T
    const x_min::T
    const k0::Int
    const Nc::Int
    const species::Vector{SpeciesState{T}}
    const a::Vector{T}
    const ex::Vector{T}
    const jy::Vector{T}
    const jz::Vector{T}
    const hprev::Vector{T}
    const chunks::Int
    const bufs::Vector{ChunkBuffers{T}}
    const kL::Int
    const kR::Int
    n::Int
    t::T
    in_left::T
    out_left::T
    out_right::T
end

function VlasovMaxwell(; Δx, Δt, N, x_min, pml::PML, source, plasma_cells::AbstractUnitRange,
                       species, scheme = PFCScheme(), chunks::Integer = 4Threads.nthreads(),
                       cut = 1e-13, probes = (pml.N + 4, N - pml.N - 2))
    T = float(promote_type(typeof(Δx), typeof(Δt)))
    op = Yee1D(; Δx = T(Δx), Δt = T(Δt), source, x_min = T(x_min), pml)
    mesh = YeeMesh1D{T}(N)
    k0, Nc = first(plasma_cells), length(plasma_cells)
    k0 ≥ pml.N + 2 && k0 + Nc ≤ N - pml.N + 1 || throw(ArgumentError(
        "plasma cells $plasma_cells reach into an absorbing layer of $(pml.N) cells " *
        "on a grid of $N"))
    chunks ≥ 1 || throw(ArgumentError("chunks = $chunks"))
    sp = SpeciesState{T}[]
    for s in species
        size(s.f, 2) == Nc || throw(DimensionMismatch(
            "species $(s.name) has $(size(s.f, 2)) space cells, the plasma $Nc"))
        push!(sp, _species_state(T, s, Nc, cut))
    end
    np = maximum(length(s.pc) for s in sp; init = 3)
    bufs = [ChunkBuffers{T}(np, Nc + 1) for _ in 1:chunks]
    kL, kR = probes
    vm = VlasovMaxwell{T, typeof(op), typeof(scheme)}(mesh, op, scheme, T(Δx), T(Δt), T(x_min),
            k0, Nc, sp, zeros(T, Nc + 1), zeros(T, Nc + 1), zeros(T, N + 1), zeros(T, N + 1),
            zeros(T, N), chunks, bufs, kL, kR, 0, zero(T), zero(T), zero(T), zero(T))
    _gauss!(vm)
    return vm
end

function _species_state(::Type{T}, s::Species, Nc, cut) where {T}
    pe = Vector{T}(s.edges)
    Δp = diff(pe)
    pc = (pe[1:end-1] .+ pe[2:end])./2
    np = length(pc)
    # The limiter coefficient of `PFCNonUniform` and its reconstruction's
    # three width sums, with the grid's own end widths standing for the empty
    # cells beyond it.
    ξ, r3, r2p, r2m = similar(Δp), similar(Δp), similar(Δp), similar(Δp)
    for j in 1:np
        d₋, d₊ = Δp[max(j - 1, 1)], Δp[min(j + 1, np)]
        ξ[j] = slope_limit(min(d₋, Δp[j], d₊)/max(d₋, Δp[j], d₊))
        r3[j] = 1/(d₋ + Δp[j] + d₊)
        r2p[j] = 1/(Δp[j] + d₊)
        r2m[j] = 1/(d₋ + Δp[j])
    end
    f = zeros(T, np, Nc + 4)
    f[:, 3:Nc+2] .= s.f
    jlo, jhi = fill(np + 1, Nc + 4), zeros(Int, Nc + 4)
    for c in 1:Nc
        nz = findall(!iszero, view(s.f, :, c))
        isempty(nz) || ((jlo[c+2], jhi[c+2]) = (first(nz), last(nz)))
    end
    fmax = maximum(s.f; init = zero(T))
    nzc = findall(c -> jlo[c+2] ≤ jhi[c+2], 1:Nc)
    lo, hi = isempty(nzc) ? (1, 0) : (first(nzc), last(nzc))
    q, m = T(s.charge), T(s.mass)
    # cells below the edge at 0 move left, the rest right; without such an
    # edge the cell that straddles 0 would move both ways, which is refused
    k0 = findfirst(iszero, pe)
    k0 === nothing && throw(ArgumentError("species $(s.name): the momentum grid needs an edge at 0"))
    return SpeciesState{T}(s.name, q, m, q/m, pe, pc, Δp, 1 ./ Δp, pe.^2, pc.^2, ξ, r3, r2p, r2m, k0,
        fmax, T(cut)*fmax, f, zeros(T, np, Nc + 4), jlo, jhi, fill(np + 1, Nc + 4),
        zeros(Int, Nc + 4), zeros(T, Nc + 1), zeros(T, 2, Nc), lo, hi, 1, 0, 1, 0,
        zero(T), zero(T), zero(T), zero(T), zero(T))
end

# Gauss's law from the left end of the plasma's grid, where `Eₓ = 0`.
function _gauss!(vm::VlasovMaxwell{T}) where {T}
    fill!(vm.ex, zero(T))
    for c in 1:vm.Nc
        ρ = zero(T)
        for s in vm.species
            ρ += s.q*_column_density(s, c)
        end
        vm.ex[c+1] = vm.ex[c] + vm.Δx*ρ
    end
    return vm
end

@inline function _column_density(s::SpeciesState{T}, c) where {T}
    n = zero(T)
    col = c + 2
    @inbounds @simd for j in s.jlo[col]:s.jhi[col]
        n += s.f[j, col]*s.Δp[j]
    end
    return n
end

# ------------------------------------------------------------------ threads

# The `k`-th of `n` contiguous chunks of `r`.
@inline function _chunk(r::UnitRange{Int}, k, n)
    len = length(r)
    q, rem = divrem(len, n)
    start = first(r) + (k - 1)*q + min(k - 1, rem)
    return start:(start + q - 1 + (k ≤ rem ? 1 : 0))
end

# `body(k, range)` over the chunks of `r`, on the threads when there is more
# than one chunk and enough work for them; serially, and without allocating,
# otherwise. The chunks are the same whichever way they run, so the answer is.
@inline function _chunked(body::F, r::UnitRange{Int}, n::Int) where {F}
    if n == 1 || length(r) < 4n
        for k in 1:min(n, max(length(r), 1))
            body(k, _chunk(r, k, min(n, max(length(r), 1))))
        end
    else
        Threads.@threads :dynamic for k in 1:n
            body(k, _chunk(r, k, n))
        end
    end
    return nothing
end

# ------------------------------------------------------------------ the step

"""
    step!(vm::VlasovMaxwell)

One step of `vm`, from `t` to `t + Δt`:

1. `a^{n+1/2} = a^{n-1/2} − Δt E_y^n` on the plasma's faces, and `(κa)²` there;
2. `X(Δt/2)`, `P(Δt)`, `X(Δt/2)` for each species, `Eₓ` advanced by the charge
   each `X` moves;
3. the transverse current, the mean of those before and after `P`, into the
   Yee step to `E_y^{n+1}`, `B_z^{n+3/2}`.
"""
function step!(vm::VlasovMaxwell{T}) where {T}
    Δt = vm.Δt
    k0 = vm.k0
    ey = vm.mesh.ey
    @inbounds for e in eachindex(vm.a)
        vm.a[e] -= Δt*ey[k0 + e - 1]
    end
    for s in vm.species
        @inbounds for e in eachindex(vm.a)
            s.A2[e] = (s.κ*vm.a[e])^2
        end
    end
    h = Δt/2
    for s in vm.species
        xsweep!(vm, s, h)
    end
    for s in vm.species
        psweep!(vm, s, Δt)
    end
    for s in vm.species
        xsweep!(vm, s, h)
    end
    _current!(vm)
    copyto!(vm.hprev, vm.mesh.hz)
    advance!(vm.mesh, vm.op, vm.t + Δt, (y = vm.jy, z = vm.jz))
    vm.n += 1
    vm.t = vm.n*Δt
    _probe!(vm)
    return vm
end

# γ − 1 without the cancellation of `√(1 + u) − 1` for small `u`.
@inline _γm1(u) = u/(1 + sqrt(1 + u))

# --------------------------------------------------- the PFC reconstruction

# The limited slopes of the PFC reconstruction towards the right neighbour `g`
# (`ϵ⁺`) and from the left one (`ϵ⁻`). As `_nu_ϵ⁺`/`_nu_ϵ⁻` in the advection
# module, the lower bound 0, but with the room to the upper bound taken as zero
# where a cell is already past it: a split sweep of a compressible flow can
# lift a cell average over the initial maximum, and the unguarded limiter then
# turned the slope round.
@inline function _ϵr(f, g, ξ, hi)
    f < g ? min(g - f, ξ*max(f, zero(f))) : max(g - f, -ξ*max(hi - f, zero(f)))
end
@inline function _ϵl(f, g, ξ, hi)
    f < g ? max(f - g, -ξ*max(f, zero(f))) : min(f - g, ξ*max(hi - f, zero(f)))
end

# The mass in the rightmost `α` of a cell (`0 ≤ α ≤ d₀`), on a non-uniform grid
# of widths `d₋, d₀, d₊`: `PFCNonUniform`'s outgoing flux through a right face,
# with the width sums' reciprocals `r3 = 1/(d₋+d₀+d₊)`, `r2p = 1/(d₀+d₊)`,
# `r2m = 1/(d₋+d₀)` given.
@inline _right_mass(α, f₋, f₀, f₊, d₋, d₀, d₊, ξ, hi, r3, r2p, r2m) =
    α*(f₀ + (d₀ - α)*r3*(_ϵr(f₀, f₊, ξ, hi)*r2p*(d₀ + d₋ - α) +
                         _ϵl(f₀, f₋, ξ, hi)*r2m*(d₊ + α)))

# The same on a uniform grid, `c` in cells, as mass over the cell width: the
# signed flux through the face between cells 0 and 1 (`c > 0`, cell 0 upwind)
# or between −1 and 0 (`c < 0`, cell 0 upwind).
@inline _flux_right(c, f₋, f₀, f₊, hi) =
    c*(f₀ + (1 - c)/3*(_ϵr(f₀, f₊, 2one(c), hi)/2*(2 - c) + _ϵl(f₀, f₋, 2one(c), hi)/2*(1 + c)))
@inline _flux_left(c, f₋, f₀, f₊, hi) =
    c*(f₀ - (1 + c)/3*(_ϵr(f₀, f₊, 2one(c), hi)/2*(1 - c) + _ϵl(f₀, f₋, 2one(c), hi)/2*(2 + c)))

# ------------------------------------------------------------------ supports

@inline _union(a1, b1, a2, b2) = (min(a1, a2), max(b1, b2))

# Shrink `lo:hi` of column `col` of `A` to its non-zero rows; an empty column
# is `(size(A, 1) + 1, 0)`, which every union takes as nothing.
@inline function _trim(A, col, lo, hi)
    @inbounds while lo ≤ hi && iszero(A[lo, col])
        lo += 1
    end
    @inbounds while hi ≥ lo && iszero(A[hi, col])
        hi -= 1
    end
    return lo ≤ hi ? (lo, hi) : (size(A, 1) + 1, 0)
end

# Zero the rows of column `col` of `A` in `olo:ohi` outside `nlo:nhi`.
@inline function _clear_outside!(A, col, olo, ohi, nlo, nhi)
    T = eltype(A)
    @inbounds for j in olo:min(ohi, nlo - 1)
        A[j, col] = zero(T)
    end
    @inbounds for j in max(olo, nhi + 1):ohi
        A[j, col] = zero(T)
    end
    return nothing
end

# The columns of `g` outside `c1:c2` that may hold values are cleared, and `f`
# and `g` change places with their supports; `c1:c2` are what was written.
function _swap!(s::SpeciesState, c1, c2)
    np = length(s.pc)
    @inbounds for c in s.glo:s.ghi
        c1 ≤ c ≤ c2 && continue
        col = c + 2
        _clear_outside!(s.g, col, s.gjlo[col], s.gjhi[col], 1, 0)
        s.gjlo[col], s.gjhi[col] = np + 1, 0
    end
    s.f, s.g = s.g, s.f
    s.jlo, s.gjlo = s.gjlo, s.jlo
    s.jhi, s.gjhi = s.gjhi, s.jhi
    s.glo, s.ghi = s.lo, s.hi
    s.lo, s.hi = c1, c2
    _shrink!(s)
    return nothing
end

# Drop empty cells from either end of `lo:hi`.
function _shrink!(s::SpeciesState)
    @inbounds while s.lo ≤ s.hi && s.jlo[s.lo + 2] > s.jhi[s.lo + 2]
        s.lo += 1
    end
    @inbounds while s.hi ≥ s.lo && s.jlo[s.hi + 2] > s.jhi[s.hi + 2]
        s.hi -= 1
    end
    return nothing
end

# ------------------------------------------------------------------ x sweep

"""
    xsweep!(vm, s, h)

Species `s` carried in `x` for a time `h`, by the PFC flux through every face
of the cells it may reach, at the velocity `(γ(p⁺) − γ(p⁻))/Δp` of each
momentum cell on that face. `Eₓ` at each face loses the charge that crossed it.
"""
xsweep!(vm::VlasovMaxwell{T, OP, PFCScheme}, s::SpeciesState, h) where {T, OP} =
    _xsweep_pfc!(vm, s, h)

function _xsweep_pfc!(vm::VlasovMaxwell{T}, s::SpeciesState{T}, h) where {T}
    s.lo > s.hi && return nothing
    Nc = vm.Nc
    c1, c2 = max(s.lo - 1, 1), min(s.hi + 1, Nc)
    for b in vm.bufs
        fill!(b.acc, zero(T))
    end
    λ = T(h)/vm.Δx
    _chunked(c1:c2, vm.chunks) do k, r
        _xchunk!(vm.bufs[k], s, vm.ex, r, c2, λ, vm.Δx, Nc)
    end
    for b in vm.bufs
        s.lost_x += b.acc[1]
        s.lost_x_energy += b.acc[2]
    end
    _swap!(s, c1, c2)
    return nothing
end

function _xchunk!(buf, s::SpeciesState{T}, ex, r, c2, λ, Δx, Nc) where {T}
    f, g, jlo, jhi = s.f, s.g, s.jlo, s.jhi
    rng = buf.rng
    Φl, Φr = buf.Φa, buf.Φb
    llo, lhi = rng[1], rng[2]
    rlo, rhi = rng[3], rng[4]
    # face `first(r)`: its rows are the supports of the two cells beside it
    e = first(r)
    _clear!(Φl, llo, lhi)
    llo, lhi = _union(jlo[e+1], jhi[e+1], jlo[e+2], jhi[e+2])
    _xflux!(Φl, buf.tmp, s, e, λ, llo, lhi)
    @inbounds for c in r
        col = c + 2
        e = c + 1
        _clear!(Φr, rlo, rhi)
        rlo, rhi = _union(jlo[e+1], jhi[e+1], jlo[e+2], jhi[e+2])
        _xflux!(Φr, buf.tmp, s, e, λ, rlo, rhi)
        nlo, nhi = _union(llo, lhi, rlo, rhi)
        nlo, nhi = _union(nlo, nhi, jlo[col], jhi[col])
        @simd for j in nlo:nhi
            g[j, col] = f[j, col] + (Φl[j] - Φr[j])
        end
        _clear_outside!(g, col, s.gjlo[col], s.gjhi[col], nlo, nhi)
        s.gjlo[col], s.gjhi[col] = _trim(g, col, nlo, nhi)
        _face!(ex, buf.acc, Φl, s, c, llo, lhi, Δx, Nc)
        Φl, Φr = Φr, Φl
        llo, lhi, rlo, rhi = rlo, rhi, llo, lhi
    end
    last(r) == c2 && _face!(ex, buf.acc, Φl, s, c2 + 1, llo, lhi, Δx, Nc)
    # leave the buffers as Φa, Φb with their rows recorded
    if Φl === buf.Φa
        rng[1], rng[2], rng[3], rng[4] = llo, lhi, rlo, rhi
    else
        rng[1], rng[2], rng[3], rng[4] = rlo, rhi, llo, lhi
    end
    return nothing
end

@inline function _clear!(Φ, lo, hi)
    @inbounds for j in lo:hi
        Φ[j] = zero(eltype(Φ))
    end
end

# The flux through face `e`, the left face of cell `e` (column `e + 2`), for the
# momentum rows `lo:hi`, in units of `f` times a cell width; at velocity
# `(γ(p_{j+1}) − γ(p_j))/Δp_j` on the face.
#
# `G` is scratch for `γ − 1` at the edges `lo:hi+1`, one square root each. The
# velocity has the sign of the momentum, `γ` being even in it and growing with
# it, so the rows below the edge at 0 (`j0`) move left and the rest right, and
# each half takes its own upwind formula without a branch.
@inline function _xflux!(Φ, G, s::SpeciesState, e, λ, lo, hi)
    lo > hi && return Φ
    f, pe2, iΔp, A2 = s.f, s.pe2, s.iΔp, s.A2[e]
    fm = s.fmax
    @inbounds @simd for i in lo:hi+1
        G[i] = _γm1(pe2[i] + A2)
    end
    j0 = clamp(s.j0, lo, hi + 1)
    @inbounds @simd for j in lo:j0-1
        c = λ*(G[j+1] - G[j])*iΔp[j]
        Φ[j] = _flux_left(c, f[j, e+1], f[j, e+2], f[j, e+3], fm)
    end
    @inbounds @simd for j in j0:hi
        c = λ*(G[j+1] - G[j])*iΔp[j]
        Φ[j] = _flux_right(c, f[j, e], f[j, e+1], f[j, e+2], fm)
    end
    return Φ
end

# Face `e`'s share of a sweep: `Eₓ` there loses the charge that crossed it, and
# at the two ends of the grid what crossed leaves the run and is counted.
@inline function _face!(ex, acc, Φ, s::SpeciesState{T}, e, lo, hi, Δx, Nc) where {T}
    Δp = s.Δp
    n = zero(T)
    @inbounds @simd for j in lo:hi
        n += Φ[j]*Δp[j]
    end
    @inbounds ex[e] -= s.q*Δx*n
    if e == 1 || e == Nc + 1
        ε = zero(T)
        A2 = s.A2[e]
        @inbounds for j in lo:hi
            ε += abs(Φ[j])*Δp[j]*_γm1(s.pc2[j] + A2)
        end
        acc[1] += abs(n)*Δx
        acc[2] += s.m*ε*Δx
    end
    return nothing
end

# ------------------------------------------------------------------ p sweep

"""
    psweep!(vm, s, Δt)

Species `s` carried in momentum for `Δt` under `F = κEₓ − (γ(x⁺) − γ(x⁻))/Δx`
at every momentum edge of every cell, `Eₓ` the mean of the cell's two faces.

Conservative semi-Lagrangian: the new average of a cell is the mass between
the feet `p − FΔt` of its two edges, the mass below a foot being the
cumulative mass up to its cell less the reconstruction's mass above it in the
cell. No Courant limit, so the finest cells may move by many widths; what
passes the ends of the grid is counted. Only the edges whose feet bracket the
column's support are visited, found by bisection, the feet being monotone in
the edge. Values below `cut` are flushed and counted. The current's sums `w`
are taken before and after.
"""
function psweep!(vm::VlasovMaxwell{T}, s::SpeciesState{T}, Δt) where {T}
    s.wlo, s.whi = s.lo, s.hi
    s.lo > s.hi && return nothing
    for b in vm.bufs
        fill!(b.acc, zero(T))
    end
    iΔx = 1/vm.Δx
    _chunked(s.lo:s.hi, vm.chunks) do k, r
        buf = vm.bufs[k]
        @inbounds for c in r
            _pcolumn!(vm.scheme, buf, s, vm.ex, c, T(Δt), iΔx, vm.Δx)
        end
    end
    for b in vm.bufs
        s.lost_p += b.acc[3]
        s.lost_p_energy += b.acc[4]
        s.flushed += b.acc[5]
    end
    _swap!(s, s.lo, s.hi)
    return nothing
end

# The foot of momentum edge `e` of cell `c` (faces `c`, `c + 1`).
@inline _foot(s, e, E, A2l, A2r, Δt, iΔx) =
    @inbounds s.pe[e] - (E - (_γm1(s.pe2[e] + A2r) - _γm1(s.pe2[e] + A2l))*iΔx)*Δt

# The edges `ea:eb` whose feet bracket `pe[jlo]` and `pe[jhi+1]`: `ea` the last
# with its foot at or below the first, or 1; `eb` the first at or above the
# second, or `np + 1`.
@inline function _bracket(s, E, A2l, A2r, Δt, iΔx, jlo, jhi)
    np = length(s.pc)
    pl, ph = s.pe[jlo], s.pe[jhi+1]
    # ea: largest e with foot(e) ≤ pl
    if _foot(s, 1, E, A2l, A2r, Δt, iΔx) > pl
        ea = 1
    else
        a, b = 1, np + 1      # foot(a) ≤ pl; find the last such
        _foot(s, b, E, A2l, A2r, Δt, iΔx) ≤ pl && (a = b)
        while b - a > 1
            m = (a + b) >>> 1
            _foot(s, m, E, A2l, A2r, Δt, iΔx) ≤ pl ? (a = m) : (b = m)
        end
        ea = a
    end
    # eb: smallest e with foot(e) ≥ ph
    if _foot(s, np + 1, E, A2l, A2r, Δt, iΔx) < ph
        eb = np + 1
    else
        a, b = ea, np + 1     # foot(b) ≥ ph; find the first such
        _foot(s, a, E, A2l, A2r, Δt, iΔx) ≥ ph && (b = a)
        while b - a > 1
            m = (a + b) >>> 1
            _foot(s, m, E, A2l, A2r, Δt, iΔx) ≥ ph ? (b = m) : (a = m)
        end
        eb = b
    end
    return ea, eb
end

# The current's two sums of a column, `Σ (f + g) Δp/γ` with `γ` at its left and
# right faces: the column before the momentum sweep, `f` with support
# `flo:fhi`, and after it, `g` with `glo:ghi`, in one pass over both.
@inline function _wsums(s::SpeciesState{T}, col, flo, fhi, glo, ghi, A2l, A2r) where {T}
    f, g = s.f, s.g
    wl = wr = zero(T)
    lo, hi = _union(flo, fhi, glo, ghi)
    @inbounds @simd for j in lo:hi
        μ = (f[j, col] + g[j, col])*s.Δp[j]
        wl += μ/(1 + _γm1(s.pc2[j] + A2l))
        wr += μ/(1 + _γm1(s.pc2[j] + A2r))
    end
    return wl, wr
end

function _pcolumn!(::PFCScheme, buf, s::SpeciesState{T}, ex, c, Δt, iΔx, Δx) where {T}
    f, g = s.f, s.g
    col = c + 2
    jlo, jhi = s.jlo[col], s.jhi[col]
    np = length(s.pc)
    A2l, A2r = s.A2[c], s.A2[c+1]
    if jlo > jhi
        s.w[1, c] = s.w[2, c] = zero(T)
        _clear_outside!(g, col, s.gjlo[col], s.gjhi[col], 1, 0)
        s.gjlo[col], s.gjhi[col] = np + 1, 0
        return nothing
    end
    pe, Δp, ξ = s.pe, s.Δp, s.ξ
    C, M = buf.C, buf.M
    hi = s.fmax
    E = s.κ*(ex[c] + ex[c+1])/2
    # C[e] the mass below edge e, for the edges of the support
    C[jlo] = zero(T)
    @inbounds for j in jlo:jhi
        C[j+1] = C[j] + f[j, col]*Δp[j]
    end
    total = C[jhi+1]
    ea, eb = _bracket(s, E, A2l, A2r, Δt, iΔx, jlo, jhi)
    pl, ph = pe[jlo], pe[jhi+1]
    kk = jlo
    @inbounds for e in ea:eb
        foot = _foot(s, e, E, A2l, A2r, Δt, iΔx)
        if foot ≤ pl
            M[e] = zero(T)
        elseif foot ≥ ph
            M[e] = total
        else
            while pe[kk+1] ≤ foot
                kk += 1
            end
            while pe[kk] > foot
                kk -= 1
            end
            f₀ = f[kk, col]
            f₋ = kk > 1 ? f[kk-1, col] : zero(T)
            f₊ = kk < np ? f[kk+1, col] : zero(T)
            d₀ = Δp[kk]
            d₋ = kk > 1 ? Δp[kk-1] : d₀
            d₊ = kk < np ? Δp[kk+1] : d₀
            M[e] = C[kk+1] - _right_mass(pe[kk+1] - foot, f₋, f₀, f₊, d₋, d₀, d₊, ξ[kk], hi,
                                         s.r3[kk], s.r2p[kk], s.r2m[kk])
        end
    end
    _pwrite!(buf, s, M, total, c, ea, eb, jlo, jhi, A2l, A2r, Δx, true)
    return nothing
end

# The new column, rows `ea:eb-1`, from the masses `M` below the feet: flushed
# below `cut` (`positive`: everything below it, else only `|f| < cut`), the
# current's sums completed, and what left through the momentum ends counted.
@inline function _pwrite!(buf, s::SpeciesState{T}, M, total, c, ea, eb, jlo, jhi, A2l, A2r,
                          Δx, positive) where {T}
    g, iΔp = s.g, s.iΔp
    col = c + 2
    cut = s.cut
    flushed = zero(T)
    @inbounds for j in ea:eb-1
        μ = M[j+1] - M[j]
        fn = μ*iΔp[j]
        if positive ? fn < cut : abs(fn) < cut
            flushed += μ
            fn = zero(T)
        end
        g[j, col] = fn
    end
    _clear_outside!(g, col, s.gjlo[col], s.gjhi[col], ea, eb - 1)
    lo, hi = _trim(g, col, ea, eb - 1)
    s.gjlo[col], s.gjhi[col] = lo, hi
    s.w[1, c], s.w[2, c] = _wsums(s, col, jlo, jhi, lo, hi, A2l, A2r)
    below, above = M[ea], total - M[eb]
    buf.acc[3] += (below + above)*Δx
    # their energy, at the grid's end cells, where f is a tail if anything
    buf.acc[4] += s.m*(below*_γm1(s.pc2[1] + A2l) + above*_γm1(s.pc2[end] + A2l))*Δx
    buf.acc[5] += flushed*Δx
    return nothing
end

# ------------------------------------------------------------------ current

# The Yee step's current argument, `−j_y Δt`, at the plasma's nodes: each
# face takes half of each neighbouring cell's `Σ f Δp/γ` towards it, averaged
# over the two sums of the step, before and after the momentum sweep.
function _current!(vm::VlasovMaxwell{T}) where {T}
    jy = vm.jy
    fill!(jy, zero(T))
    k0, Δt, a = vm.k0, vm.Δt, vm.a
    for s in vm.species
        s.wlo > s.whi && continue
        coef = Δt*s.q*s.κ/4
        @inbounds for c in s.wlo:s.whi
            jy[k0 + c - 1] += coef*a[c]*s.w[1, c]
            jy[k0 + c] += coef*a[c+1]*s.w[2, c]
        end
    end
    return nothing
end

# ------------------------------------------------------------------ fluxes

# The incoming and outgoing electromagnetic energy, integrated in time at the
# two probe nodes: `F²` and `G²` at the left one, `F² − G² = E_y B_z` at the
# right one, `B_z` at the node and the step's end the mean of its two cells
# at the two half steps either side.
function _probe!(vm::VlasovMaxwell{T}) where {T}
    ey, hz, hp = vm.mesh.ey, vm.mesh.hz, vm.hprev
    k = vm.kL
    E = ey[k]
    B = (hz[k-1] + hz[k] + hp[k-1] + hp[k])/4
    F, G = (E + B)/2, (E - B)/2
    vm.in_left += F*F*vm.Δt
    vm.out_left += G*G*vm.Δt
    k = vm.kR
    E = ey[k]
    B = (hz[k-1] + hz[k] + hp[k-1] + hp[k])/4
    vm.out_right += E*B*vm.Δt
    return nothing
end

# ------------------------------------------------------------------ spline

# Moments `y2[k1:k2]` of the clamped (zero end slope) cubic spline through
# `(x[k], y[k])`, `k1:k2`, with spacings `h[k] = x[k+1] − x[k]`; `d` is scratch
# for the pivots. A Thomas sweep.
function _spline_moments!(y2, d, y, h, k1, k2)
    @inbounds begin
        # row k1: 2h m + h m'; interior: h₋ m₋ + 2(h₋ + h) m + h m₊; row k2: h₋ m₋ + 2h₋ m
        hk = h[k1]
        d[k1] = 2hk
        y2[k1] = 6*((y[k1+1] - y[k1])/hk)
        for k in k1+1:k2
            h₋ = h[k-1]
            diag = k < k2 ? 2(h₋ + h[k]) : 2h₋
            rhs = k < k2 ? 6*((y[k+1] - y[k])/h[k] - (y[k] - y[k-1])/h₋) :
                           -6*((y[k] - y[k-1])/h₋)
            μ = h₋/d[k-1]
            d[k] = diag - μ*h₋
            y2[k] = rhs - μ*y2[k-1]
        end
        y2[k2] /= d[k2]
        for k in k2-1:-1:k1
            y2[k] = (y2[k] - h[k]*y2[k+1])/d[k]
        end
    end
    return y2
end

# The spline at `x` in the interval `[xk, xk + hk]`.
@inline function _spline_at(x, xk, hk, yk, yk1, mk, mk1)
    A = (xk + hk - x)/hk
    B = 1 - A
    return A*yk + B*yk1 + ((A*A*A - A)*mk + (B*B*B - B)*mk1)*(hk*hk)/6
end

# The spline sweep in momentum: the primitive over the support and eight cells
# either side, where it is flat, splined and differenced between the feet.
function _pcolumn!(::SplineScheme, buf, s::SpeciesState{T}, ex, c, Δt, iΔx, Δx) where {T}
    f, g = s.f, s.g
    col = c + 2
    jlo, jhi = s.jlo[col], s.jhi[col]
    np = length(s.pc)
    A2l, A2r = s.A2[c], s.A2[c+1]
    if jlo > jhi
        s.w[1, c] = s.w[2, c] = zero(T)
        _clear_outside!(g, col, s.gjlo[col], s.gjhi[col], 1, 0)
        s.gjlo[col], s.gjhi[col] = np + 1, 0
        return nothing
    end
    pe, Δp = s.pe, s.Δp
    C, M, y2, d = buf.C, buf.M, buf.y2, buf.tmp
    E = s.κ*(ex[c] + ex[c+1])/2
    k1, k2 = max(1, jlo - 8), min(np + 1, jhi + 9)
    @inbounds for k in k1:jlo
        C[k] = zero(T)
    end
    @inbounds for j in jlo:jhi
        C[j+1] = C[j] + f[j, col]*Δp[j]
    end
    total = C[jhi+1]
    @inbounds for k in jhi+2:k2
        C[k] = total
    end
    _spline_moments!(y2, d, C, Δp, k1, k2)
    ea, eb = _bracket(s, E, A2l, A2r, Δt, iΔx, jlo, jhi)
    kk = k1
    pl, ph = pe[k1], pe[k2]
    @inbounds for e in ea:eb
        foot = _foot(s, e, E, A2l, A2r, Δt, iΔx)
        if foot ≤ pl
            M[e] = zero(T)
        elseif foot ≥ ph
            M[e] = total
        else
            while pe[kk+1] ≤ foot
                kk += 1
            end
            while pe[kk] > foot
                kk -= 1
            end
            M[e] = _spline_at(foot, pe[kk], Δp[kk], C[kk], C[kk+1], y2[kk], y2[kk+1])
        end
    end
    _pwrite!(buf, s, M, total, c, ea, eb, jlo, jhi, A2l, A2r, Δx, false)
    return nothing
end

xsweep!(vm::VlasovMaxwell{T, OP, SplineScheme}, s::SpeciesState, h) where {T, OP} =
    _xsweep_spline!(vm, s, h)

# The spline sweep in `x`, row by row: each momentum row of the cells the
# species may reach, with two empty cells either side, is gathered, its
# primitive splined, and differenced between the feet of the faces.
function _xsweep_spline!(vm::VlasovMaxwell{T}, s::SpeciesState{T}, h) where {T}
    s.lo > s.hi && return nothing
    Nc = vm.Nc
    np = length(s.pc)
    c1, c2 = max(s.lo - 1, 1), min(s.hi + 1, Nc)
    g1, g2 = c1 - 2, c2 + 2           # the gathered cells; ghosts are empty
    ncell = g2 - g1 + 1
    nknot = ncell + 1
    f, g, Δp = s.f, s.g, s.Δp
    Δx = vm.Δx
    # the rows any of them holds
    rlo, rhi = np + 1, 0
    for c in c1-1:c2+1
        rlo, rhi = _union(rlo, rhi, s.jlo[c+2], s.jhi[c+2])
    end
    for c in c1:c2
        col = c + 2
        _clear_outside!(g, col, s.gjlo[col], s.gjhi[col], rlo, rhi)
    end
    for b in vm.bufs
        fill!(b.acc, zero(T))
        fill!(view(b.jx, 1:Nc+1), zero(T))
        fill!(view(b.xpiv, 1:nknot), Δx)
    end
    _chunked(rlo:rhi, vm.chunks) do k, r
        buf = vm.bufs[k]
        C, M, y2, d, hx = buf.C, buf.M, buf.y2, buf.tmp, buf.xpiv
        @inbounds for j in r
            C[1] = zero(T)
            for i in 1:ncell
                C[i+1] = C[i] + f[j, g1 + i + 1]*Δx      # the column of cell g1+i-1
            end
            _spline_moments!(y2, d, C, hx, 1, nknot)
            for e in c1:c2+1
                A2 = s.A2[e]
                v = (_γm1(s.pe2[j+1] + A2) - _γm1(s.pe2[j] + A2))*s.iΔp[j]
                xk = (e - g1)*Δx                 # face e relative to knot 1
                foot = xk - v*h
                kk = clamp(floor(Int, foot/Δx) + 1, 1, nknot - 1)
                Mf = _spline_at(foot, (kk - 1)*Δx, Δx, C[kk], C[kk+1], y2[kk], y2[kk+1])
                M[e - c1 + 1] = C[e - g1 + 1] - Mf      # mass through face e, rightwards
            end
            for c in c1:c2
                g[j, c + 2] = f[j, c + 2] + (M[c - c1 + 1] - M[c - c1 + 2])/Δx
            end
            for e in c1:c2+1
                Φ = M[e - c1 + 1]
                buf.jx[e] += Φ*Δp[j]
                if e == 1 || e == Nc + 1
                    buf.acc[1] += abs(Φ)*Δp[j]
                    buf.acc[2] += s.m*abs(Φ)*Δp[j]*_γm1(s.pc2[j] + s.A2[e])
                end
            end
        end
    end
    for b in vm.bufs
        @inbounds for e in c1:c2+1
            vm.ex[e] -= s.q*b.jx[e]
        end
        s.lost_x += b.acc[1]
        s.lost_x_energy += b.acc[2]
    end
    for c in c1:c2
        col = c + 2
        s.gjlo[col], s.gjhi[col] = _trim(g, col, rlo, rhi)
    end
    _swap!(s, c1, c2)
    return nothing
end

# ------------------------------------------------------------------ diagnostics

"""
    node_x(vm), cell_x(vm)

The coordinates of the Yee nodes and of the plasma's cell centres.
"""
node_x(vm::VlasovMaxwell) = vm.x_min .+ vm.Δx .* (0:vm.mesh.N)
cell_x(vm::VlasovMaxwell) = vm.x_min .+ vm.Δx .* ((vm.k0 - 1) .+ (0.5:1:vm.Nc))

"""
    plasma_nodes(vm)

The Yee nodes that are the faces of the plasma's cells, `k0 : k0 + Nc`.
"""
plasma_nodes(vm::VlasovMaxwell) = vm.k0:(vm.k0 + vm.Nc)

"""
    distribution(vm, s)

The distribution of species `s` (an index), `[momentum cell, plasma cell]`, a
view. Zero outside each column's support, `support(vm, s, c)`.
"""
distribution(vm::VlasovMaxwell, s::Integer) = view(vm.species[s].f, :, 3:vm.Nc+2)

"""
    support(vm, s, c)

The momentum cells of plasma cell `c` outside which species `s` is zero.
"""
support(vm::VlasovMaxwell, s::Integer, c::Integer) =
    vm.species[s].jlo[c+2]:vm.species[s].jhi[c+2]

"""
    density!(n, vm, s), density(vm, s)

The density of species `s` over the plasma's cells, `Σ f Δp`, in `n_c`.
"""
function density!(n, vm::VlasovMaxwell, s::Integer)
    sp = vm.species[s]
    for c in 1:vm.Nc
        n[c] = _column_density(sp, c)
    end
    return n
end
density(vm::VlasovMaxwell{T}, s::Integer) where {T} = density!(zeros(T, vm.Nc), vm, s)

"""
    kinetic_energy(vm[, s])

`Σ m f (γ − 1) Δp Δx` of species `s`, or of all, in `n_c mₑc² (c/ω₀)`, with
`γ` the mean of the two faces' at the step's end, `a^n` the mean of the two
half steps' (`a^{n+1/2} = a^{n-1/2} − Δt E_y^n` estimated from the current
field).
"""
function kinetic_energy(vm::VlasovMaxwell{T}, si::Integer) where {T}
    s = vm.species[si]
    s.lo > s.hi && return zero(T)
    K = zero(T)
    k0 = vm.k0
    ey = vm.mesh.ey
    for c in s.lo:s.hi
        al = vm.a[c] - vm.Δt/2*ey[k0 + c - 1]
        ar = vm.a[c+1] - vm.Δt/2*ey[k0 + c]
        Al, Ar = (s.κ*al)^2, (s.κ*ar)^2
        col = c + 2
        @inbounds for j in s.jlo[col]:s.jhi[col]
            p2 = s.pc2[j]
            K += s.f[j, col]*s.Δp[j]*(_γm1(p2 + Al) + _γm1(p2 + Ar))/2
        end
    end
    return s.m*K*vm.Δx
end
kinetic_energy(vm::VlasovMaxwell{T}) where {T} =
    sum(kinetic_energy(vm, s) for s in eachindex(vm.species); init = zero(T))

"""
    longitudinal_energy(vm)

`Σ Eₓ²/2 Δx` over the plasma's faces.
"""
longitudinal_energy(vm::VlasovMaxwell) = sum(abs2, vm.ex)*vm.Δx/2

"""
    field_energy(vm[, k1, k2])

The transverse field's energy between nodes `k1` and `k2` (by default the two
probes): `Σ E_y²/2 Δx` over the nodes and `Σ B_z^{n-1/2} B_z^{n+1/2}/2 Δx` over
the cells between them, the energy Yee's update conserves.
"""
function field_energy(vm::VlasovMaxwell{T}, k1 = vm.kL, k2 = vm.kR) where {T}
    ey, hz, hp = vm.mesh.ey, vm.mesh.hz, vm.hprev
    U = zero(T)
    for k in k1:k2
        U += ey[k]^2
    end
    for k in k1:k2-1
        U += hz[k]*hp[k]
    end
    return U*vm.Δx/2
end

"""
    energy_account(vm)

The energy budget between the two probes, a NamedTuple in `n_c mₑc² (c/ω₀)`:
`field` (transverse, [`field_energy`](@ref)), `longitudinal`, `kinetic`,
`lost` (particles that left through the ends of `x` or `p`, at their energy),
`out_left` (`∫G²dt` at the left probe), `out_right` (`∫E_yB_z dt` at the right
one), `in_left` (`∫F²dt` at the left probe), and `residual`, the sum of the
first six less `in_left`, which a conservative run keeps small.
"""
function energy_account(vm::VlasovMaxwell{T}) where {T}
    field = field_energy(vm)
    longitudinal = longitudinal_energy(vm)
    kinetic = kinetic_energy(vm)
    lost = sum(s.lost_x_energy + s.lost_p_energy for s in vm.species; init = zero(T))
    total = field + longitudinal + kinetic + lost + vm.out_left + vm.out_right
    return (; field, longitudinal, kinetic, lost, out_left = vm.out_left,
            out_right = vm.out_right, in_left = vm.in_left, residual = total - vm.in_left)
end

end # module
