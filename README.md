# percMPI
Percolation for node-distributed data

How to run
==========
Following packages are required
    CMake
    MPI (testing on MPICH2)

Visual Studio Code
==================
Settings for Visual Studio Code are included (.vscode)
Runs out-of-the-box on a linux machine with the following packages installed:
    visual studio code
    gcc
    gdb
    CMake
    MPI
    clang-format

Using shell-commands
==================
When started from the main folder
scripts/build_make.sh builds and compiles different versions
parallel or multiple: Version for multiple nodes
parallel-collectives or multiple-collectives: Version for multiple nodes with MPI Collectives where possible
sequential or single: Version for a single node, the sequential algorithm

