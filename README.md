# Vasilek — Vlasov Adaptive SImulator of pLasma Electrodynamics and Kinetics

An ongoing project on developing of a parallel 2D2P Maxwell — Vlasov — Boltzmann solver on adaptive meshes.

As for now, the following funcionality has been impemented:
* PFC scheme on 1D non-uniform grid
* Strang splitting for 1D1V simulations
* 1D Poisson Fourier solver

The program has been verified on the following tests:
* [Long-lasting plasma oscillations](verification/plasma-oscillations-1d1v.ipynb)
* [Landau damping](verification/landau-damping-1d1v.ipynb)