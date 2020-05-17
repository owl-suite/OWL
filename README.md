# OWL
Open-source/Oak-Ridge Wang-Landau (OWL): A suite for first-principles based Monte Carlo simulations

OWL is a scientific software for large-scale Monte Carlo simulations for materials. Originally developed at Oak Ridge National Laboratory to implement Wang-Landau sampling with first-principles calculations for the study of finite temperature materials properties, it is now open-sourced and has extended to provide a collection of commonly used classical, modern and parallel Monte Carlo algorithms.

OWL is written in C++ with an object-oriented, modular software architecture, adopting the "MPI+X" programming model for parallelization. It provides two modes of simulation: the stand-alone mode for simulating user-implemented model Hamiltonians; and the driver mode that utilizes an external package as a library for the calculations of physical observables. Today, OWL interfaces with two open-source density functional theory codes, Quantum Espresso and Locally Self-Consistent Multiple Scattering, to perform first-principles based statistical mechanics simulations. It is therefore highly suitable for running (flying) on high performance computers.

OWL is under active development with proper software engineering practices. Capability extensions are driven by ongoing research activities. Our goal is to build a community code and create a platform to facilitate algorithm advancements and knowledge exchanges. We welcome comments, suggestions and contributions from the scientific community.
