// Created by Huzaifa Mustafa Unjhawala

#ifndef FE_H
#define FE_H

// All the required libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <Eigen>
//#include <Eigen/Dense>



typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
class FE{
    public:
        FE(unsigned int nelx,unsigned int nely,unsigned int nelz,double length, double breadth,double height,double youngs_mod, double pois_rat,double density); // The constructor takes the required arguments
        double xi_at_node(unsigned int local_node_no); // Used in basis function and basis gradient of 27 node element
        double eta_at_node(unsigned int local_node_no);// Used in basis function and basis gradient of 27 node element
        double kap_at_node(unsigned int local_node_no);// Used in basis function and basis gradient of 27 node element
        double basis_function(unsigned int node, double xi,double eta, double kap); //Calculates the basis function corresponding to the node "node" and at the point 'xi','eta' and 'kap' in the parametric space
        std::vector<double> basis_gradient(unsigned int node, double xi,double eta,double kap); //Calculates the gradient of the basis function - similar to above
        void mesh(unsigned int no_quad_points); // Function to mesh the domain - Fills the Nodal connectivity and the Elemental Conncectivity matrices - Can even handle different number of elements along each axis
        void define_boundary_condition(double g1,double g2); // Function to fill the boundary_values (stores the values at the boundaries) and the boundary_nodes (stores the global node number of the nodes on the boundary) - The force is not used since there is no force to apply for this problem
        double C(unsigned int i, unsigned int j, unsigned int k, unsigned int l); // Used to get the elasticity tensor C
        void init_data_structs(); //To resize all the global matrices based on the mesh - Internal function
        void cal_jac(unsigned int q1, unsigned int q2,unsigned int q3);
        void cal_k_local(); // Calculates the K local for one element - As all elements are the same, can use the same klocal
        void cal_m_local();// Calculates the K local for one element - As all elements are the same, can use the same mlocal
        void assemble(); //Uses the klocal to assemble K global using the volume fractions
        Eigen::VectorXd solve(); // Solves and returns U which is then used in the toplogy code
        void initial_cond(); // Apply the initial conditions
        void solve_trans(double gamma,double beta,unsigned int time_steps,double delta_t); // SOlve for transient cases
        void fem_to_vtk(); // Write to vtk for steady case
        void write_to_vtk(int t_step); // write to vtk for dynamic case


        // Class datastructures
        double L,B,H,g1_,g2_,f1,E,nu,lambda,mu,penal_,detJ,rho; //Standard constants - L - Length, B - breadth, g1 - Dirichlet conditon, E - Youngs Modulus, nu - Poissons ration, lambda and mu
        // are the lame's parameters
        unsigned int nelx_,nely_,nelz_,nel,nnx_,nny_,nnz_,no_of_nodes,no_of_nodes_per_element,total_dofs,dofs_per_ele,quad_rule,dim; // standard mesh descriptions


        std::vector<std::vector<double> > NC; //Nodal Connectivity - NC[i] gives the x,y and z - coordinate of the ith global node. Size - (No.of nodes, dim)
    // klocal and mlocal matrices defined as eigen matrices for better speed
        Eigen::MatrixXd Klocal;
        Eigen::MatrixXd Mlocal;
        std::vector<std::vector<unsigned int> > EC_2; //Elemental connectivity - EC[i][j] gives the global node number for local node j in element i - Size - (No. of elements, No. of nodes per element)
        std::vector<double> boundary_values; // Vector having the dirichlet boundary value wherever its defined and 0 in all other entries - Size (No. of nodes)
        std::vector< unsigned int > boundary_nodes; //Vector having all the nodes that are part of the dirichlet boundary - Size depends on the number of dirichlet nodes
        std::vector<std::vector<double> > quad_points; // Vector for the location of the quad points in terms of xi
        std::vector<double> quad_weights; // Vector for the weights at the quad points. Still 1D as the weigths do not depend on xi or eta
        Eigen::MatrixXd invJ; // Inverse Jacobian needed
        Eigen::VectorXd U; // Using Eigne vector to define to solution for ease of solving the linear equation
        Eigen::VectorXd F; // No forcing but the dirichlet conditions will apply
        Eigen::MatrixXd K; // Global stiffness matrix
        Eigen::MatrixXd M; // Global mass matrix
        Eigen::MatrixXd damp; // Global C matrix (damping)
        Eigen::VectorXd u_n; // This is the eigen vector to define the previous d in our time loop
        Eigen::VectorXd v_n;  // The velocity vector
        Eigen::VectorXd a_n; // The acceleration vector
};

#endif

