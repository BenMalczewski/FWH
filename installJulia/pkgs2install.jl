#=
    _______       ____  __  |
   / ____/ |     / / / / /  |
  / /_   | | /| / / /_/ /   |   Ffowcs-Williams Hawking Solver
 / __/   | |/ |/ / __  /    |   Version: 1.0.0
/_/      |__/|__/_/ /_/     |   Written By: Benjamin Malczewski & Sam Salehian, Ph.D.

Script to install all of the required packages for FWH
=#

using Pkg


# Install Required Packages
Pkg.add("MPI")
Pkg.add("MPIPreferences")
Pkg.add("PencilArrays")
Pkg.add("FFTW")
Pkg.add("Memoization")
Pkg.add("LinearAlgebra")
Pkg.add("Statistics")
Pkg.add("DelimitedFiles")

