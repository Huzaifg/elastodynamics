#!/usr/bin/env zsh

path_to_eigen="$PWD/Eigen"
if [[ -d "$res" ]]
then
	echo "Compiling code for the parameters specified in main.cpp"
	g++ -I $path_to_eigen main.cpp FE.cpp -Wall -O3 -o run
	echo "Compilation complete! Now running executable"
	# Run with nelx,nely,nelz,length,breadth,height,dirichlet boundary condition,gamma,beta,delta T,number of time steps
	./run 2 2 20 0.1 0.1 1 0.05 0.5 0.25 0.025 200
else
	mkdir res
	echo "Compiling code for the parameters specified in main.cpp"
	g++ -I $path_to_eigen main.cpp FE.cpp -Wall -O3 -o run
	echo "Compilation complete! Now running executable"
	# Run with nelx,nely,nelz,length,breadth,height,dirichlet boundary condition,gamma,beta,delta T,number of time steps
	./run 2 2 20 0.1 0.1 1 0.05 0.5 0.25 0.025 200
fi