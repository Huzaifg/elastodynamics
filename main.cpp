#include "FE.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>

#define NSEC_PER_SEC 1000000000



int main(int argc, char* argv[]){
    struct timespec start, end;
    if(argc < 11){
        std::cout<<"Please provide sufficient inputs to run the elastodynamics code - You need 11 inputs"<<std::endl<<
        "nelx,nely,nelz,length,breadth,height,dirichlet boundary condition,gamma,beta,delta T,number of time steps"<<std::endl;
        return 0;
    }
    char *pCh;
	unsigned int no_quad_points = 3; // Do not change
    unsigned int nelx = strtoul(argv[1], &pCh, 10); // Number of elements along x
    unsigned int nely = strtoul(argv[2], &pCh, 10); // Number of elements along y
    unsigned int nelz = strtoul(argv[3], &pCh, 10);
    double length = std::stod(argv[4]); // Length
    double breadth = std::stod(argv[5]); // Breadth
    double height = std::stod(argv[6]);
//    double youngs_mod = 2 * pow(10,11);
    double youngs_mod = 1000; // Can be changed for different materials but before compiling
    double pois_rat = 0.3; // Can be changed for different materials but before compiling
    double g1 = 0.; // Can be changed for different materials but before compiling
    double g2 = std::stod(argv[7]);
    double gamma = std::stod(argv[8]);
    double beta = std::stod(argv[9]);
    double delta_t = std::stod(argv[10]);
    unsigned int time_steps = strtoul(argv[11], &pCh, 10);
    double density = 1.; //Can be changed for different materials but before compiling
    FE tr(nelx,nely,nelz,length,breadth,height,youngs_mod,pois_rat,density);
    tr.mesh(no_quad_points);
    
    
    tr.init_data_structs();
    tr.define_boundary_condition(g1,g2);
    tr.cal_k_local();
    tr.cal_m_local();
    clock_gettime(CLOCK_MONOTONIC, &start);
    tr.assemble();
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_t elapsed_sec_1 = end.tv_sec - start.tv_sec;
    long elapsed_nsec_1 = end.tv_nsec - start.tv_nsec;

    double elapsed_total_1 = elapsed_sec_1 + (double)elapsed_nsec_1 / (double)NSEC_PER_SEC;
//
    printf("Time taken for fe_imp %g \n",elapsed_total_1 * 1000);
    Eigen::VectorXd U;
    clock_gettime(CLOCK_MONOTONIC, &start);
    U = tr.solve();
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_t elapsed_sec_2 = end.tv_sec - start.tv_sec;
    long elapsed_nsec_2 = end.tv_nsec - start.tv_nsec;
//
    double elapsed_total_2 = elapsed_sec_2 + (double)elapsed_nsec_2 / (double)NSEC_PER_SEC;
//
    printf("Time taken for solution of linear system %g \n",elapsed_total_2 * 1000);
    
    tr.fem_to_vtk(); // Writes the steady state to a vtk
    tr.solve_trans(gamma,beta,time_steps,delta_t); // soves for transient, writes to vtk within the function itself
}
