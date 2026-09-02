# Normalization conventions

Every quantity in this package is dimensionless. The conventions were, until
now, recorded only in the prose of the verification notebooks — and the
`PoissonFourier1D` bug that returned a field too large by `1/Δx²` came directly
from that gap: the module and the notebooks disagreed about what a wavenumber
meant, and nothing said which was right.

## Base units

For the 1D1V electrostatic problem, with `Nₑ` the equilibrium electron density,
`m` and `e` the electron mass and elementary charge:

| quantity | normalized to | symbol |
|---|---|---|
| time | inverse plasma frequency | `ωₚ⁻¹`, where `ωₚ² = 4πe²Nₑ/m` |
| velocity | thermal velocity | `v_th` |
| length | `v_th/ωₚ` | the Debye length |
| density | `Nₑ` | |
| electric field | `m ωₚ v_th / e` | |
| distribution function | `Nₑ/v_th` | |

With those, the Vlasov–Poisson system is

```
∂f/∂t + v ∂f/∂x + E ∂f/∂v = 0
∂E/∂x = nᵢ - nₑ,       nₑ = ∫ f dv
```

and the total energy the verification runs report is

```
ε = ∬ f v²/2 dv dx + ∫ E²/2 dx
```

For the electromagnetic case (`wakefield.jl`) time and length are normalized to
the laser frequency and `c/ω₀` instead, and the field to `m c ω₀ / e`, so that
the dimensionless laser amplitude is `a₀`.

## Wavenumbers

**`PoissonFourier1D` takes a physical wavenumber**, `k = 2πm/L` with `L = NΔx`,
not a per-sample one. This is the distinction that produced the bug:
`rfftfreq(n)` returns cycles per *sample*, so `Δx` has to be supplied as the
sampling rate — `rfftfreq(n, 1/Δx)` — or it never enters the spectrum at all
while the φ→E difference divides by it once.

The scheme satisfies an exact discrete identity on a grid mode, which the test
asserts to machine precision rather than to a tolerance:

```
E_num[i] == (sin(kΔx)/(kΔx)) · E_exact[x_i]
```

The `sinc` factor is the centred difference; the DFT contributes nothing
further, being exact on grid modes.

## The advection argument

`advect!(dest, src, scheme, c, ws)` takes, as its fourth argument:

* the **Courant number** `vΔt/Δx` for the uniform-grid schemes;
* the **displacement** `vΔt`, a length, for `PFCNonUniform`.

A non-uniform grid has no single Courant number to quote, so there is no
common dimensionless form. This is a wart; it is documented rather than hidden.

## Currents in the FDTD solver

`make_advance_fields` adds its current argument **straight into the field**:

```julia
f.ey[i] += jy[i]
```

so the caller owes it the time step. The argument is `-J Δt`, not `J`. Getting
this wrong is not subtle in its consequences — the wakefield example reached a
peak field of `1.0e22` against a source amplitude of `1.0` before it was
corrected.

Sign convention, taken from the longitudinal direction that the
plasma-oscillation notebook verifies against the analytic plasma frequency:

```
∂p/∂t = +E        ∂E/∂t = -n·p
```

Together these give an oscillation. Flipping either gives exponential growth.

The frequency of that oscillation is now asserted rather than asserted-about:
`test/test_verification.jl` measures it at `ω = 1.005719` against the
Bohm–Gross `√(1 + 3k²) = 1.005904` at `k = 2π/100`, which is 0.018% — and
0.57% from the cold `ωₚ = 1`, so the thermal correction is resolved rather than
merely tolerated. Until that test was written this paragraph was the only place
the claim existed.

## Boundaries in the FDTD solver

Both ends of the grid are **perfect electric conductors**. The solver imposes
that by never writing to either end node: the interior loop runs
`pml.N+2 : Nx-pml.N` and the two absorbing-layer loops stop short of both ends,
so `ey[1]`, `ez[1]`, `ey[end]` and `ez[end]` hold the zero `YeeMesh1D` gives
them, for all time.

They are not dead storage. `update_hz!` reads `ey[end]`, which is how the
condition enters the solution, and a pulse reaching either wall comes back with
its sign inverted — measured reflection coefficient −0.9998. It is also what
makes the staggered energy exactly conserved: the discrete curls are adjoint
only because the boundary terms vanish.

Two consequences for callers:

* **Seed the interior.** Writing a nonzero value into an end node does not
  launch a wave, it silently changes the boundary condition to `E = const`.
* **The current arrays span all `N+1` nodes, but only `2:N` are driven.** A
  current cannot be injected into a perfect conductor. The first and last
  entries of `j.y` and `j.z` are ignored; they exist so the arrays index
  alongside `f.ey`.

The second of those was wrong until recently: the current was added over
`1:Nx`, so node 1 was driven while node `Nx+1` was not. It broke the symmetry
the energy identity rests on, and injected a source into a conductor at one end
only.

## Bounds for the PFC schemes

`PFC` and `PFCNonUniform` require `fmin` and `fmax`. By Liouville's theorem `f`
stays within the extrema of the **initial** condition along characteristics, so
these are global initial bounds, supplied by the caller. They cannot be derived
per line inside a multidimensional sweep, which is why there is no default —
one silently wrong default already cost a factor of 740 between two overloads
that were supposed to agree.
