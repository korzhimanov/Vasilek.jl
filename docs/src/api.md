```@meta
CurrentModule = Vasilek
```

# API reference

Every docstring in the package: the exported names, and the internals that the
docstrings and the notes refer to. The build fails if one is missing here.

## Workspaces

One generic function serves every kind of object that needs scratch memory:
schemes, the Strang step, the Poisson solver and the collision operators.

```@docs
workspace
```

## Advection

```@docs
Advection
AbstractAdvection1D
advect!
```

### Schemes

```@docs
Upwind
LaxWendroff
Godunov
SemiLagrangian
PFC
PFCNonUniform
```

### Options

```@docs
Advection.AbstractReconstruction
PiecewiseConstant
PiecewiseLinear
Advection.AbstractLimiter
NoLimiter
VanLeer
Superbee
Advection.AbstractSpline
LinearSpline
QuadraticSpline
CubicSpline
```

### Internals

```@docs
Advection.slope_limit
Advection._validate
Advection._validate_workspace
Advection._check_bounds
Advection._validate_courant
```

## Strang splitting

```@docs
StrangSplitting.strang_step!
StrangSplitting.StrangWorkspace
StrangSplitting.make_time_step_2d!
```

## The 1D1V Vlasov–Poisson driver

```@docs
VlasovPoisson1D1V
vlasov_poisson
```

### Internals

```@docs
VlasovPoisson1D1V.make_poisson
VlasovPoisson1D1V.cell_widths
VlasovPoisson1D1V.nlogn
VlasovPoisson1D1V.line_advector
VlasovPoisson1D1V.substeps
VlasovPoisson1D1V.mode_amplitude
VlasovPoisson1D1V.keeps_bounds
```

## Poisson solver

```@docs
PoissonFourier1D.PoissonFFT1D
PoissonFourier1D.solve!
PoissonFourier1D.generate_solver
```

## Maxwell solver

```@docs
FDTD1D.YeeMesh1D
FDTD1D.PML
FDTD1D.make_advance_fields
```

## Collisions

```@docs
Collisions
AbstractCollisionOperator
collide!
BGK
Collisions.Landau1P
```

### Internals

```@docs
Collisions._discrete_maxwellian!
Collisions._maxwellian_moments
Collisions._sampled_maxwellian!
Collisions.∂f∂v
```
