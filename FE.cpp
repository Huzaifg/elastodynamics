// Created by Huzaifa Mustafa Unjhawala 

#include "FE.h"
#include <fstream>
#include <iostream>




// Constructor

FE::FE(unsigned int nelx, unsigned int nely,unsigned int nelz, double length, double breadth,double height, double youngs_mod, double pois_rat,double density){
	L = length;
	B = breadth;
    H = height;
	nelx_ = nelx;
	nely_ = nely;
    nelz_ = nelz;
	E = youngs_mod;
	nu = pois_rat;
    rho = density;
}
// Basis function - Internal function needed for fe implementation
inline double FE::basis_function(unsigned int node , double xi, double eta, double kap){
//Hardcoded for the 8 node.
    double output;
    switch(node) {
        case 0:
            output = 0.125 * (1-xi) * (1-eta) * (1-kap);
            break;
        case 1:
            output = 0.125 * (1+xi) * (1-eta) * (1-kap);
            break;
        case 2:
            output = 0.125 * (1+xi) * (1+eta) * (1-kap);
            break;
        case 3:
            output = 0.125 * (1-xi) * (1+eta) * (1-kap);
            break;
        case 4:
            output = 0.125 * (1-xi) * (1-eta) * (1+kap);
            break;
        case 5:
            output = 0.125 * (1+xi) * (1-eta) * (1+kap);
            break;
        case 6:
            output = 0.125 * (1+xi) * (1+eta) * (1+kap);
            break;
        case 7:
            output = 0.125 * (1-xi) * (1+eta) * (1+kap);
            break;
        default:
            std::cout<<"There is no "<<node<<" node, invalid input"<<std::endl;
            break;
    }
    return output;
}
// Basis function defined using the general formula obtained from
inline std::vector<double> FE::basis_gradient(unsigned int node,double xi, double eta,double kap){
// Kind off hard coded graident for 8 node
    
    std::vector<double> bg(dim,0.0);
    switch(node) {
        case 0:
            bg[0] = -0.125 * (1 - eta) * (1-kap);
            bg[1] = -0.125 * (1 - xi) * (1-kap);
            bg[2] = -0.125 * (1 - xi) * (1-eta);
            break;
        case 1:
            bg[0] = 0.125 * (1 - eta) * (1-kap);
            bg[1] = -0.125 * (1 + xi) * (1-kap);
            bg[2] = -0.125 * (1 + xi) * (1-eta);
            break;
        case 2:
            bg[0] = 0.125 * (1 + eta) * (1 - kap);
            bg[1] = 0.125 * (1 + xi) * (1-kap);
            bg[2] = -0.125 * (1 + xi) * (1+eta);
            break;
        case 3:
            bg[0] = -0.125 * (1 + eta) * (1-kap);
            bg[1] = 0.125 * (1 - xi) * (1-kap);
            bg[2] = -0.125 * (1 - xi) * (1+eta);
            break;
        case 4:
            bg[0] = -0.125 * (1 - eta) * (1+kap);
            bg[1] = -0.125 * (1 - xi) * (1+kap);
            bg[2] = 0.125 * (1 - xi) * (1-eta);
            break;
        case 5:
            bg[0] = 0.125 * (1 - eta) * (1+kap);
            bg[1] = -0.125 * (1 + xi) * (1+kap);
            bg[2] = 0.125 * (1 + xi) * (1-eta);
            break;
        case 6:
            bg[0] = 0.125 * (1 + eta) * (1+kap);
            bg[1] = 0.125 * (1 + xi) * (1+kap);
            bg[2] = 0.125 * (1 + xi) * (1+eta);
            break;
        case 7:
            bg[0] = -0.125 * (1 + eta) * (1+kap);
            bg[1] = 0.125 * (1 - xi) * (1+kap);
            bg[2] = 0.125 * (1 - xi) * (1+eta);
            break;
        default:
            std::cout<<"There is no "<<node<<" node, invalid input"<<std::endl;
    }
    return bg;
}
// Function to generate mesh
void FE::mesh(unsigned int no_quad_points){
	std::cout<<"Generating Mesh .."<<std::endl;

    // The number of nodes is just an extension of 2D
    no_of_nodes = (nelx_ + 1) * (nely_ + 1) * (nelz_ + 1);
    std::cout<<"Total no. of nodes is "<<no_of_nodes<<std::endl;
    // Each node has 3 degrees of freedom as the number of dimensions is 3
    dim = 3;
    total_dofs = no_of_nodes * dim;

    // Nodal Coordinate array will give the coordinates of each dof not just each node
    NC.resize(total_dofs);
    // Since NC is a vector of a vector , we need to initialize each row of EC
    for(int i = 0; i < total_dofs; i++){
        NC[i] = std::vector<double>(dim,0.0);
    }



    // Make NC

    double incr_x = L/(nelx_); // since the nodes are equally spaced, the x coordinate will differ by this increment
    double incr_y = B/(nely_); // similarly, the y coordinate will differ by this incremenet
    double incr_z = H/(nelz_);
    double x = 0.0; // first node is 0,0
    double y = 0.0;
    double z = 0.0;
    // Construct NC - NC[i][0] gives the x - coordinate of the ith global node, NC[i][1] gives the y and NC[i][2] gives the z
    // Here, 3 dofs make up one node and hence pairs of dofs will have the same coordinates
    for(int i = 0; i < total_dofs - 1 ; i = i + dim){
        NC[i][0] = x;
        NC[i+1][0]  = x;
        NC[i+2][0]  = x;
        x += incr_x;
        NC[i][1] = y;
        NC[i+1][1] = y;
        NC[i+2][1] = y;
        NC[i][2] = z;
        NC[i+1][2] = z;
        NC[i+2][2] = z;
        // if we reach the x and y limit, increment z and set x and y to 0
        if((abs(NC[i][0] - L) < 0.0000001) && (abs(NC[i][1] - B) < 0.0000001)){
            x = 0;
            y = 0;
            z += incr_z;
        }
        // If we have reached the x limit, reset x  to 0 and increment y
        else if(abs(NC[i][0] - L) < 0.0000001){
            x = 0;
            y += incr_y;
        }
    }
    no_of_nodes_per_element = 8;
    // Since each node has more than dim dofs, the dofs per element will be the no of nodes per element * dimensions
    dofs_per_ele = no_of_nodes_per_element * dim;
    nel = nelx_ * nely_ * nelz_;
    //Element connectivity matrix
    EC_2.resize(nel);
    for(int i = 0; i < nel; i++){
        // Over here we have to use dofs_per_ele as these will be the number of columns in EC
        EC_2[i] = std::vector<unsigned int>(no_of_nodes_per_element);
    }
    
    
    int nnx_ = nelx_ + 1; // Number of nodes along x
    int nny_ = nely_ + 1; // Number of nodes along y
    int nnz_ = nelz_ + 1; // Number of nodes along z
    int inc_x_ = 0; // Tells how many increments we have had in x
    int inc_y_ = 0; // Tells how many increments we have had in y
    int inc_z_ = 0; // Tells how many increments we have has in z
    
    // Construct EC - EC[i][j] gives the global node number for local node j in element i
    for(int i = 0; i < nel;i++){
        //If we have reached last node on x, we increment y and reset our x coutner
        if(inc_x_ == nnx_ - 1){
            inc_x_ = 0;
            inc_y_ += 1;
            // if this increment causes us to reach the last node on y as well then we need to increment z
            if(inc_y_ == nny_ - 1){
                inc_z_ += (nnx_ + nny_ - 1); // incrememnt z by the number of nodes we have passed already
                inc_y_ = 0;
            }
            
        }

//            Storing clockwise on each element - pattern was hand derived using some examples for lesser number of elements
//            Node numbers increase left to right. Inc_y takes us to the nodes of the element 1 level up and inc_z takes us to the nodes of element 1 level deeper
        EC_2[i][0] = i + inc_y_ + inc_z_ ;
        EC_2[i][1] = i + 1 + inc_y_ + inc_z_ ;
        EC_2[i][2] = nnx_ + 1 + i + inc_y_ + inc_z_ ;
        EC_2[i][3] = nnx_ + i + inc_y_ + inc_z_ ;
        EC_2[i][4] = i + (nnx_ * nny_) + inc_y_ + inc_z_ ;
        EC_2[i][5] = i + (nnx_ * nny_) + 1 + inc_y_ + inc_z_ ;
        EC_2[i][6] = i + (nnx_ * nny_) + 1 + nnx_ + inc_y_ + inc_z_;
        EC_2[i][7] = i + (nnx_ * nny_)  + nnx_ + inc_y_ + inc_z_ ;
        inc_x_ += 1;
    }
    // Set up quadrature data - Cant change number of quad points for now - Can include functionality with simple if-else if needed
    quad_rule = no_quad_points;
    quad_points.resize(quad_rule); //Resize quadpoints to appropriate size
    for(int i = 0; i < quad_rule; i++){
        quad_points[i] = std::vector<double>(dim);
    }
    quad_weights.resize(quad_rule); //Resize quadweights to appropriate size

    quad_points[0][0] = -sqrt(3./5.); // xi
    quad_points[0][1] = -sqrt(3./5.); // eta
    quad_points[0][2] = -sqrt(3./5.); // kappa
    quad_points[1][0] = 0.; // xi
    quad_points[1][1] = 0.; // eta
    quad_points[1][2] = 0.; // kappa
    quad_points[2][0] = sqrt(3./5.); // xi
    quad_points[2][1] = sqrt(3./5.); // eta
    quad_points[2][2] = sqrt(3./5.); // kappa

    quad_weights[0] = 5./9. ;
    quad_weights[1] = 8./9. ;
    quad_weights[2] = 5./9. ;

}

void FE::define_boundary_condition(double g1,double g2){
    std::cout<<std::endl<<"Defining boundary condition"<<std::endl;
//    Initialize the Dirichlet and Neumann boundary condtions
    g1_ = g1;
    g2_ = g2;
//    At each dof which is a boundary, we will put the value of g1, else we will put 0
    boundary_values.resize(total_dofs);
    for(unsigned int dof_no = 0; dof_no < total_dofs ; dof_no = dof_no + 3){
        // Checks the z coordinate of a patricular dof to see if its 0
        if(NC[dof_no][2] == 0){
            boundary_values[dof_no] = g1_;
            boundary_nodes.push_back(dof_no);
            boundary_values[dof_no+1] = g1_;
            boundary_nodes.push_back(dof_no+1);
            boundary_values[dof_no+2] = g1_;
            boundary_nodes.push_back(dof_no+2);
        }
        // Checks the z coordinate of a patricular dof to see if its H
        else if(abs(NC[dof_no][2] - H) <0.00001){
//            boundary_values[dof_no] = 0.; // Only Y displacement
//            boundary_nodes.push_back(dof_no);
            boundary_values[dof_no+1] = g2_; // Only Y displacement
            boundary_nodes.push_back(dof_no+1);
//            boundary_values[dof_no+2] = 0.; // Only Y displacement
//            boundary_nodes.push_back(dof_no+2);
        }
    }
}

// Function used for debigging - exports the eigen matrices as csv files for easy visualization
void saveData(std::string fileName, Eigen::MatrixXd  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
 
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

// Function to initilize all our data structures
void FE::init_data_structs(){
    std::cout<<"Initializing data structures"<<std::endl;
    K.resize(total_dofs,total_dofs); //Resize K
    K.setZero(total_dofs,total_dofs); // Initialize K to 0
    M.resize(total_dofs,total_dofs);
    M.setZero();
    damp.resize(total_dofs,total_dofs);
    damp.setZero();
    F.resize(total_dofs); //Resize F
    F.setZero(total_dofs); // Setting F to zero here itself since we know the size
    U.resize(total_dofs); //Resive d
    u_n.resize(total_dofs);
    v_n.resize(total_dofs);
    a_n.resize(total_dofs);
}

// Function for calculating the value of C - elasticity tensor
inline double FE::C(unsigned int i, unsigned int j, unsigned int k, unsigned int l){
    double lambda = (E * nu)/((1. + nu) * (1. - 2.*nu));
    double mu = E/(2. *(1. + nu));
    return lambda * (i==j) * (k==l) + mu * ((i==k)*(j==l) + (i==l)* (j==k));
}

// This function is used to evaluate the Jacobian matrix, its inverse and its determinant at a paticular guassian quadrature point.
inline void FE::cal_jac(unsigned int q1, unsigned int q2,unsigned int q3){
    Eigen::MatrixXd Jac;
    Jac.resize(dim,dim);
    invJ.resize(dim,dim);
    for(unsigned int i = 0; i < dim; i++){
        for(unsigned int j = 0; j < dim; j++){
            Jac(i,j) = 0;
            // Looping through the nodes of an element
            for(unsigned int A = 0; A < no_of_nodes_per_element; A ++){
                // Over here dim*A is used because EC has dim dofs per node. Each of these dofs have the same coordinate, so we can pick either one while calculating the jacobian. Over here, we use all the even dofs
                Jac(i,j) += NC[dim*EC_2[0][A]][i] * basis_gradient(A, quad_points[q1][0], quad_points[q2][1],quad_points[q3][2])[j];
            }
        }
    }
    detJ = Jac.determinant();
    invJ = Jac.inverse();
    saveData("invJ.csv", invJ);
}


//This method is used to fill up the elemental K matrix.Since our domain is made up of all the same type of elements with equal sizing, the cal_k_local method is only called once to find the K elemental for the 1st element since all the elements will have the same k local.
void FE::cal_k_local(){
    std::cout<<"Determining Klocal"<<std::endl;
    // Initializing Klocal
    Klocal.resize(dofs_per_ele,dofs_per_ele);
    Klocal.setZero();
//    for(int res = 0; res < dofs_per_ele; res++){
//        Klocal[res] = std::vector<double>(dofs_per_ele);
//    }
    

//    std::fill(Klocal.begin(), Klocal.end(), std::vector<double>(dofs_per_ele, 0.));
    for(unsigned int q1 = 0; q1 < quad_rule ; q1++){
        for(unsigned int q2 = 0; q2 < quad_rule ; q2++){
            for(unsigned int q3 = 0; q3 < quad_rule; q3++){
                cal_jac(q1,q2,q3);
                // Now we go ahead and fill in the Klocal array
                for(unsigned int A = 0; A < no_of_nodes_per_element; A++){
                    // Capital I and K denote the physical coordinates
                    for(unsigned int I = 0; I < dim; I++){
                        for(unsigned int B=0 ; B < no_of_nodes_per_element; B++){
                            for(unsigned int K = 0; K < dim; K++){
                                for(unsigned int J = 0; J < dim; J++){
                                    for(unsigned int L = 0; L < dim; L++){
                                        // Looping over the parametric coordinates - I think we only need to loop over j and k since only those indicies are used - Not sure though
                                        for(unsigned int j = 0; j < dim; j++){
                                            for(unsigned int l = 0; l < dim; l++){
                                                // Added i and k since we maybe do need it - Need to figure out how to reduce these number of loops - Will be too slow
                                                Klocal(dim*A + I,dim*B + K) += (basis_gradient(A, quad_points[q1][0], quad_points[q2][1],quad_points[q2][2])[j] * invJ(j,J)) * C(I,J,K,L) * (basis_gradient(B, quad_points[q1][0], quad_points[q2][1],quad_points[q2][2])[l] * invJ(l,L)) * detJ * quad_weights[q1] * quad_weights[q2] * quad_weights[q3];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                    }
                }
            }
        }
    }
    saveData("k_local.csv", Klocal);
}
//This method is used to fill up the elemental M matrix. Since our domain is made up of all the same type of elements with equal sizing, the cal_m_local method is only called once to find the M elemental for the 1st element since all the elements will have the same m local.
void FE::cal_m_local(){
    std::cout<<"Determining M local"<<std::endl;
    Mlocal.resize(dofs_per_ele,dofs_per_ele);
    Mlocal.setZero();
    for(unsigned int q1 = 0; q1 < quad_rule ; q1++){
        for(unsigned int q2 = 0; q2 < quad_rule ; q2++){
            for(unsigned int q3 = 0; q3 < quad_rule; q3++){
                for(unsigned int A = 0; A < no_of_nodes_per_element ; A++){
                    for(unsigned int I = 0; I < dim ; I++){
                        
                        for(unsigned int B = 0; B < no_of_nodes_per_element ; B++){
                            for(unsigned int J = 0; J<dim ; J++){
    //                        Here we use quad_points[q1][0] instead of quad_points[q1][i] as we dont have an i loop. Even for K loop we can do this, it will not make a difference.
                                Mlocal(dim*A + I,dim*B + J) += rho*basis_function(A, quad_points[q1][0], quad_points[q2][1],quad_points[q3][2])*basis_function(B, quad_points[q1][0], quad_points[q2][1],quad_points[q3][2])*abs(detJ)*quad_weights[q1]*quad_weights[q2]*quad_weights[q3] * (I==J);
                            }
                        }
                    }
                }
            }
        }
    }
    saveData("m_local.csv", Mlocal);
}
// This function assembles our K amd M matrics
void FE::assemble(){
    unsigned int row;
    unsigned int col;
    std::cout<<"Assembling K and M"<<std::endl;
    for(int ele = 0; ele < nel ; ele++){
        // Now we assemble the Klocal into the K matrix which is the global matrix and similarly the mlocal
        for(unsigned int I = 0; I < no_of_nodes_per_element ; I++){
            row =dim*EC_2[ele][I];
            for(unsigned int J = 0; J < no_of_nodes_per_element ; J++){
                col =dim*EC_2[ele][J];
//                K.coeffRef(EC[ele][I],EC[ele][J]) += * Klocal(I,J);
                // Here we use heavy loop unrolling. Really speeds up the assembly part as the compiler seems a lot more code
                K(row,col) +=   Klocal(dim*I,dim*J);
                K(row,col+1) +=   Klocal(dim*I,dim*J+1);
                K(row,col+2) +=   Klocal(dim*I,dim*J+2);
                K(row+1,col) +=   Klocal(dim*I+1,dim*J);
                K(row+1,col+1) +=   Klocal(dim*I+1,dim*J+1);
                K(row+1,col+2) +=   Klocal(dim*I+1,dim*J+2);
                K(row+2,col) +=  Klocal(dim*I+2,dim*J);
                K(row+2,col+1) +=  Klocal(dim*I+2,dim*J+1);
                K(row+2,col+2) +=  Klocal(dim*I+2,dim*J+2);
                
                M(row,col) +=   Mlocal(dim*I,dim*J);
                M(row,col+1) +=   Mlocal(dim*I,dim*J+1);
                M(row,col+2) +=   Mlocal(dim*I,dim*J+2);
                M(row+1,col) +=   Mlocal(dim*I+1,dim*J);
                M(row+1,col+1) +=   Mlocal(dim*I+1,dim*J+1);
                M(row+1,col+2) +=   Mlocal(dim*I+1,dim*J+2);
                M(row+2,col) +=  Mlocal(dim*I+2,dim*J);
                M(row+2,col+1) +=  Mlocal(dim*I+2,dim*J+1);
                M(row+2,col+2) +=  Mlocal(dim*I+2,dim*J+2);
            }
        }
    }
//    std::cout << "The determinant of M is " << M.determinant() << std::endl;
//    saveData("k_before.csv", K);
//    saveData("m_before.csv", M);
    double a_damp = 1.0;
    double b_damp = 0.001;
    // Using regilegh damping to find the damping matrix C (damp)
    damp = (a_damp * M) + (b_damp * K);
    // Now we apply the Dirichlet boundary conditons and modify K accordingly
    std::cout<<"Applying Dirichlet BC's"<<std::endl;
    for(unsigned int i : boundary_nodes){
        double g = boundary_values[i];
        // Loop to move the approprate column of K to the RHS - source - https://www.math.colostate.edu/~bangerth/videos.676.21.65.html
        for(unsigned int row = 0; row < total_dofs; row++){
            // This condition is so that a dof which has already been set in F is not changed
            if(row == i){
                continue;
            }
            // All the other dofs in F are varied as we move the column of K to the RHS
            else{
                F[row] = F[row] - g * K.coeffRef(row,i);
            }
        }
        // Set all the diagonal elements to 1 and all the other elements in the row and column to 1
        K.row(i) *= 0;
        K.col(i) *= 0;
        M.row(i) *= 0;
        M.col(i) *= 0;
        damp.row(i) *= 0;
        damp.col(i) *= 0;
        K(i,i) = 1.;
        M(i,i) = 1.;
        damp(i,i) = 1.;
        // Set the value in F at athe node
        F[i] = g;
    }
//    saveData("k_after.csv", K);
//    saveData("m_after.csv", M);
//    saveData("damp_after.csv", damp);
    
}

// The solve method for the steady case. Converts K and F into sparse matrix and then solves. Very fast if converted to sparse matrix
Eigen::VectorXd FE::solve(){
    std::cout<<"Solving..."<<std::endl;
    Eigen::SparseMatrix<double> K_ = K.sparseView();
    Eigen::SparseVector<double> F_ = F.sparseView();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(K_);
    U = solver.solve(F_);
    std::cout<<U<<std::endl;
    return U;
}

// Function to apply the initial condition
void FE::initial_cond(){
//    Setting all values to 0
    u_n.setZero();
    v_n.setZero();
    
//    Applying the dirichlet values in the initial conditions
    for(unsigned int i : boundary_nodes){
        double g = boundary_values[i];
        u_n[i] = g;
    }
    
    Eigen::VectorXd RHS;
    RHS.resize(total_dofs);
    RHS = F - (damp*v_n) - (K*u_n);
    a_n = M.inverse() * RHS;

    
}

// Function to run the loop over time
void FE::solve_trans(double gamma,double beta,unsigned int time_steps,double delta_t){
    std::cout<<"Solving Transient"<<std::endl;
    // Apply inital condition
    initial_cond();
    const double delt = delta_t;
    // Initalize u tilda and v tilda
    Eigen::VectorXd u_tilde;
    u_tilde.resize(total_dofs);
    Eigen::VectorXd v_tilde;
    v_tilde.resize(total_dofs);
    

    Eigen::VectorXd RHS;
    RHS.resize(total_dofs);
    Eigen::MatrixXd inverting_m; // This is the matrix to be inverted or LHS
    // Again , we eventually convert to sparse while inverting
    Eigen::SparseMatrix<double> invm_ ;
    Eigen::SparseVector<double> RHS_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    inverting_m.resize(total_dofs,total_dofs);
    for(int t_step=0;t_step < time_steps+1; t_step++){
//        if(t_step%5 == 0){
        // Write a vtk file every time step since its pretty fast anyways
        write_to_vtk(t_step);
        // Print out the displacement field so that the user can know what is happenin g
        std::cout<<"Displacement Field "<<t_step<<std::endl;
        std::cout<<u_n<<std::endl;
//        }
        // Calculate u tilda according to formula
        u_tilde = (u_n + (delt * v_n) + (0.5 * pow(delt,2) * (1 - 2*beta)) * a_n);\
        // similarly v tilda
        v_tilde = v_n + (delt * (1 - gamma) * a_n);
        // get the rhs
        RHS = F - (damp*v_tilde) - (K*u_tilde);
        // get the matrix that needs to be inverted
        inverting_m = M + (gamma * delt * damp) + (beta * pow(delt,2) * K);
        //convert to sparse
        invm_ = inverting_m.sparseView();
        RHS_ = RHS.sparseView();
        solver.compute(invm_);
        //solve
        a_n = solver.solve(RHS_);
        
        // update u and v
        u_n = u_tilde + (0.5 * pow(delt,2) * 2 * beta * a_n);
        v_n = v_tilde + (gamma * delt * a_n);
    }
}
// Function to write to vtk for transient case - again, transcribed from matlab
void FE::write_to_vtk(int t_step){
    // Write to file all the stuff that is needed for plotting
    std::cout<<"Writing to vtu file for Transient state"<<std::endl;
    std::ofstream out_file;
    std::stringstream ss;
    ss<<"res/output_"<<t_step<<".vtu";
    std::string f_name;
    f_name+= ss.str();
//        out_file.open("res/output_" + std::to_string(scheme) + std::to_string(t_step) +  ".vtu");
    out_file.open(f_name);
    if(out_file.fail()){
        std::cout<<"File did not open"<<std::endl;
    }
//   Writing headers
    out_file<<"<?xml version=\"1.0\"?>"<<std::endl;
    out_file<<"<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >"<<std::endl;
    out_file<<"<UnstructuredGrid>"<<std::endl;
//    Writing nodal coordinates array
    out_file<<"<Piece  NumberOfPoints=\""<<no_of_nodes<<"\" NumberOfCells=\""<<nel<<"\">"<<std::endl;
    out_file<<"<Points>"<<std::endl;
    out_file<<"<DataArray type=\"Float32\" NumberOfComponents=\""<<dim<<"\" format=\"ascii\">"<<std::endl;
    float z = 0.;
    for(int node = 0; node < total_dofs; node = node+3){
        out_file<<(float) NC[node][0]<<" "<<(float) NC[node][1]<<" "<<(float) NC[node][2]<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</Points>"<<std::endl;
    
//    Writing Element connectivity
    out_file<<"<Cells>"<<std::endl;
    out_file<<"<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">"<<std::endl;
    for(int ele = 0; ele < nel; ele++){
        out_file<<EC_2[ele][0]<<" "<<EC_2[ele][1]<<" "<<EC_2[ele][2]<<" "<<EC_2[ele][3]<<" "<<EC_2[ele][4]<<" "<<EC_2[ele][5]<<" "<<EC_2[ele][6]<<" "<<EC_2[ele][7]<<std::endl;
    }

    out_file<<"</DataArray>"<<std::endl;
//    Writing element offsets vector(required in VTK format)
    unsigned int offsets = 0;
    out_file<<"<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">"<<std::endl;

    for(int ele = 0; ele < nel; ele++){
        offsets=offsets+8;
        out_file<<offsets<<std::endl;
    }


//    Writing element types vector(required in VTK format)
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">"<<std::endl;
    unsigned int ty = 12;
    for(int ele = 0; ele < nel; ele++){
        out_file<<ty<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</Cells>"<<std::endl;
//    Writing field values (if any) array
    
    out_file<<"<PointData  Scalars=\"u\">"<<std::endl;
    out_file<<"<DataArray  type=\"Float32\"  Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">"<<std::endl;
    float z1 = 0;
    for(int node = 0; node < total_dofs; node=node+3){
        out_file<<(float)u_n(node)<<" "<<(float)u_n(node+1)<<" "<<(float)u_n(node+2)<<std::endl;
//        out_file<<(float)d_n(node)<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</PointData> "<<std::endl;
    out_file<<"</Piece> "<<std::endl;
    out_file<<"</UnstructuredGrid> "<<std::endl;
    out_file<<"</VTKFile> "<<std::endl;
    out_file.close();
        
}

// function to write to vtk for steady case
void FE::fem_to_vtk(){
        
    // Write to file all the stuff that is needed for plotting
    std::cout<<"Writing to vtu file for steady state"<<std::endl;
    std::ofstream out_file;
    std::stringstream ss;
    ss<<"res/steady.vtu";
    std::string f_name;
    f_name+= ss.str();
//        out_file.open("res/output_" + std::to_string(scheme) + std::to_string(t_step) +  ".vtu");
    out_file.open(f_name);
    if(out_file.fail()){
        std::cout<<"File did not open"<<std::endl;
    }
//   Writing headers
    out_file<<"<?xml version=\"1.0\"?>"<<std::endl;
    out_file<<"<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >"<<std::endl;
    out_file<<"<UnstructuredGrid>"<<std::endl;
//    Writing nodal coordinates array
    out_file<<"<Piece  NumberOfPoints=\""<<no_of_nodes<<"\" NumberOfCells=\""<<nel<<"\">"<<std::endl;
    out_file<<"<Points>"<<std::endl;
    out_file<<"<DataArray type=\"Float32\" NumberOfComponents=\""<<dim<<"\" format=\"ascii\">"<<std::endl;
    float z = 0.;
    for(int node = 0; node < total_dofs; node = node+3){
        out_file<<(float) NC[node][0]<<" "<<(float) NC[node][1]<<" "<<(float) NC[node][2]<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</Points>"<<std::endl;
    
//    Writing Element connectivity
    out_file<<"<Cells>"<<std::endl;
    out_file<<"<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">"<<std::endl;
    for(int ele = 0; ele < nel; ele++){
        out_file<<EC_2[ele][0]<<" "<<EC_2[ele][1]<<" "<<EC_2[ele][2]<<" "<<EC_2[ele][3]<<" "<<EC_2[ele][4]<<" "<<EC_2[ele][5]<<" "<<EC_2[ele][6]<<" "<<EC_2[ele][7]<<std::endl;
    }

    out_file<<"</DataArray>"<<std::endl;
//    Writing element offsets vector(required in VTK format)
    unsigned int offsets = 0;
    out_file<<"<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">"<<std::endl;

    for(int ele = 0; ele < nel; ele++){
        offsets=offsets+8;
        out_file<<offsets<<std::endl;
    }


//    Writing element types vector(required in VTK format)
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">"<<std::endl;
    unsigned int ty = 12;
    for(int ele = 0; ele < nel; ele++){
        out_file<<ty<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</Cells>"<<std::endl;
//    Writing field values (if any) array
    
    out_file<<"<PointData  Scalars=\"u\">"<<std::endl;
    out_file<<"<DataArray  type=\"Float32\"  Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">"<<std::endl;
    float z1 = 0;
    for(int node = 0; node < total_dofs; node=node+3){
        out_file<<(float)U(node)<<" "<<(float)U(node+1)<<" "<<(float)U(node+2)<<std::endl;
//        out_file<<(float)d_n(node)<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</PointData> "<<std::endl;
    out_file<<"</Piece> "<<std::endl;
    out_file<<"</UnstructuredGrid> "<<std::endl;
    out_file<<"</VTKFile> "<<std::endl;
    out_file.close();
        
}
