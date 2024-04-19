#=
    _______       ____  __  |
   / ____/ |     / / / / /  |
  / /_   | | /| / / /_/ /   |   Ffowcs-Williams Hawking Solver
 / __/   | |/ |/ / __  /    |   Version: 1.0.0
/_/      |__/|__/_/ /_/     |   Written By: Benjamin Malczewski & Sam Salehian, Ph.D.

Script to install mpiexecjl for convenience when running FWH.
=#


using MPI


# Install mpiexecjl
MPI.install_mpiexecjl()
