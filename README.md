# Elastodynamics
This is work done for the M E 601 class - Introduction to Finite Elements class offered at the University Of Wisconsin Madison. Consider this as a novice attempt to writing FEM code. Effort was made to ensure it is as accessible as possible, however, for any clarifications, please create an issue.

## Important Information
The code solves a 3D elastostatics and elastodynamics problem for a Linear Isotropic Cantilivered Beam. Although the material porperties Youngs Modulus and Poissons Ratio are not user inputs, they can be changed in the main.cpp file prior to compiling. The code currently only supports 2 Dirichlet boundary conditions at Z=0 and Z=L. However, the value of the Dirichlet boundary condition on the z=1 face is a user input and can be varied during run time. The dimensions of the cantilivered beam and the and the number of elements along each direction can also be given during run time to analyse for different shapes.
The code supports the Newark Family of methods for time stepping and thus the user can supply the desired value of beta and gamma during run time. The user can also specify the desired Delta T and the number of time steps during run time. Be sure to cross check the stability condition while supplying a delta T.
The code also currently only supports the Bi-linear 8 node Hexahedron element. The quadrature used is a 9 point quadrature. Changes to either of these parameters will break the program.

## What does the program output?
The program first outputs a vtu file with the steady state output (steady.vtu). This file can be found in the bash res subdirectory. The code also outputs a vtu file at each time step. This can also be found in the bash res subdirectory with naming convention output_TimeStepNo.vtu. These vtu files can be visualised using ParaView or Visit.

## How to build and run?
Once the repository has been cloned and cd-ed into, there are two ways that the code can be built and run.

###1) Using CMake
```bash
mkdir build
cd build
cmake ..
make
# nelx,nelx,nelz,length,breadth,height,dirichlet boundary condition at Z=1,gamma,beta,delta T,Number of timesteps
./run 2 2 20 0.1 0.1 1 0.05 0.5 0.25 0.025 200
```
1) Using the provided shell scripts
The input paramaters can be varied in the shell script
```bash
zsh automate_it.sh
```
