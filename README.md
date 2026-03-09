# Entanglement Dynamics in 3D Anderson Model

We analyze the single-particle entanglement entropy of the 3D Anderson Model with periodic boundary conditions. In `buildAnd.c`, the Hamiltonian is diagonalized in `C` using the `SLEPc` library, and the eigenvectors and eigenvalues are calculated. The entanglement analysis and other measures are calculated in the `Julia` programming language. The aim is to give a single-parameter formulation of the entanglement dynamics; see the [paper](https://www.mdpi.com/1099-4300/28/1/29).
