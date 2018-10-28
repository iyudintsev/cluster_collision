# Clusters Collision Project

This service helps sovling the problem of Molecular Dynamics. Suppose we have a cub area with periodic boundary conditions. 
This area contains two clusters of particles.
Two clusters move to each other and a collision occurs. We choose Lennard Jones Potential to implement interaction between particles and use Verlet integration to simulate the dynamics of the system.

Besides, we use the method of neighbors search to optimize the calculations.

## How to run the simulation

Firstly, we need to download and extract Eigen library it in the project folder.

`wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2`

`tar -xf 3.3.4.tar.bz2`

`mv eigen-eigen-5a0156e40feb/ Eigen && rm 3.3.4.tar.bz2`

After that to build a project using `make`. The result of this is executable file `model`.

Run a file `viewer.py` in the folder `/veaw` to visualize the results. The dependencies are libraries `numpy` and `matplotlib`.

`python3 viewer.py`

## How to run the tests

We need to build test file `make -f MakefileTests`. As a result of it we can run executable file `test`.
