# percMPI
Percolation for node-distributed data

How to run
==========
The following packages are required
* CMake
* MPI (tested on MPICH2)

Visual Studio Code
==================
Settings for Visual Studio Code are included (.vscode).
Runs out-of-the-box on a linux machine with the following packages installed:
* visual studio code
* gcc
* gdb
* CMake
* MPI
* clang-format
    
Settings
==================
For examples on how to build and run the code, see the visual studio tasks (.vscode/tasks.json)

Set the compile flag `SINGLENODE` for running sequentially

Set the compile flag `COMMUNICATION` for running parallel

Using shell-commands
==================
When started from the main folder, `scripts/build_make.sh` builds and compiles different versions:
* `parallel` or `multiple`: Version for multiple nodes
* `parallel-collectives` or `multiple-collectives`: Version for multiple nodes with MPI Collectives where possible
* `sequential` or `single`: Version for a single node, the sequential algorithm

Data
==================
A raw data set can be loaded with `inputMode 3` (see task _isotroph_). An average and rms value need to be specified.

Random data is sampled for `inputMode 20`
