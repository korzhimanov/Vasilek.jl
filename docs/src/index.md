# Vasilek — Vlasov Adaptive Simulator of pLasma Electrodynamics and Kinetics

An ongoing project on developing a parallel 2D2P Maxwell — Vlasov — Boltzmann
solver on adaptive meshes. What exists today is 1D1V on static grids, uniform
in x and uniform or non-uniform in v; 2D2P and adaptive meshes are the goal, not
yet the code.

As for now, the following functionality has been implemented:
* Advection schemes as dispatchable types: upwind, Lax—Wendroff, Godunov
  (piecewise-constant or -linear, with flux limiters), semi-Lagrangian
  (linear, quadratic or cubic B-splines), and PFC on uniform and non-uniform grids
* Strang splitting for 1D1V simulations, and `vlasov_poisson`, the 1D1V
  electrostatic driver (optionally with a collision operator)
* 1D Poisson Fourier solver
* 1D FDTD Maxwell solver with PML
* BGK collision operator

## In this manual

* [Normalization conventions](normalization.md): the units, the sign of the
  Poisson equation, and what the fourth argument of `advect!` means.
* [Notes on the 1D1V driver](driver-notes.md): why `vlasov_poisson` makes each
  of its choices, with the measurements behind them.
* [Migrating to 0.2](migration-0.2.md): from the 0.1 closures to scheme values.
* [API](api.md): every docstring, by module, the internals included.
* [Verification](verification.md): the table of what the test suite asserts
  about the verification studies, and each study rendered with its output.
