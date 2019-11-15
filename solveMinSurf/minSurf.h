#include<iostream>
#include<stdio.h>
#include<vector>
#include<valarray>
#include<cmath>
#include<fstream>
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/Core"
#include"Eigen/Eigen/SparseCholesky"


#define EPS 1e-12


template <typename T>
T residual( const std::vector<T> &z, std::vector<T> &r ) {
    // computes residual entries in r[i]
    // returns sum of squares of r
    T sum = 0;
    for (auto& n : r)
        sum += n;
}

// ################################################################################################
// Boundary nodes by lexicographical ordering
// To be extended for a) different boundaries
//                    b) MPI?
void setBdryNodes( std::vector<int> &bdryNodeList, const int N ) {
    int i=0;
    
    // Set first and last row
    for(i=0; i<N; i++) {
    // First row
        bdryNodeList.push_back(i);
    // Last row
        bdryNodeList.push_back(N*(N-1) + i);
//         std::cout << "me: " << i << std::endl;
    }
    
    // outer edges
    for(i=N; i<N*(N-1); i+=N) {
        bdryNodeList.push_back(i);
        bdryNodeList.push_back(i+N-1);
    }

    // Debug output
    for(auto& i: bdryNodeList)
        std::cout << i << std::endl;
}

// Inner nodes by lexicographical ordering
void setInnerNodes( std::vector<int> &innerNodeList, const int N ) {
    
    // Set inner nodes (assuming structured cartesian grid)
    for(int i=1; i<N-1; i++) {
        for(int j=1; j<N-1; j++) {
            innerNodeList.push_back(i+j*N);
        }
    }
    // Debug output
    std::cout << "innernodes:\n";
    for(auto& it: innerNodeList)
        std::cout << it << std::endl;
        
}

// Function to apply boundary conditions - needs to be extended -> Sankar
template <typename mType>
void applyBC( const int BCtype, Eigen::MatrixBase<mType> &V , const std::vector<int> &bdryNodeList, 
            std::vector<int> innerNodeList ) {
    // Set boundary values
    for(auto& i: bdryNodeList)
        V[i]= sin( (double)i );
    // Set rest to zero
    for(auto& i: innerNodeList)
        V[i] = 0;
}

// ################################################################################################
template <typename mType, typename dType>
void buildPoissonMatrix( Eigen::SparseMatrix<dType> &A, const std::vector<int> &bdryNodeList, 
                         const std::vector<int> &innerNodeList, const int N ) { 

    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(3*N);

    // Assemble the FD matrix (should work, validated with matlab)
    // Set inner nodes - 5-pt stencil of FD (there should not be an issue with reaching 
    // the limits of the matrix...)
    for(auto& i: innerNodeList) {
        tripletList.push_back(triplet(i ,i    , 4));
        tripletList.push_back(triplet(i ,i+1  ,-1));
        tripletList.push_back(triplet(i ,i-1  ,-1));
        tripletList.push_back(triplet(i ,i+N-1,-1));
        tripletList.push_back(triplet(i ,i-N+1,-1));
    }

//     double v_ij=0;
//     
//     for(int i=1; i<(N-1)*(N-1); i++) {
//         for(int j=1; j<(N-1)*(N-1); j++) {
//             if (i==j)
//                 v_ij = 4;
//             else if (i-1==j and i-1>=0)
//                 v_ij = -1;
//             else if (i+1==j and i+1<N*N)
//                 v_ij = -1;
//             else if (i+N-1==j and i+N-1<N*N)
//                 v_ij = -1;
//             else if (i-N+1==j and i-N+1>=0)
//                 v_ij = -1;
//             else v_ij = 0;
//             
// //             std::cout << i << "," << j << ": " << v_ij << std::endl;
//             if (v_ij)
//                 tripletList.push_back(triplet(i,j,v_ij));
//         }
//     }

// Set a 1 where Dirichlet BC applies
    for(auto& i: bdryNodeList)
        tripletList.push_back(triplet(i, i,1));

    // Build sparse A from triplets
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout << A << std::endl;
}

template <typename mType, typename dType>
void getInitGuess( Eigen::MatrixBase<mType> &zE, const Eigen::MatrixBase<mType> &bE, const std::vector<int> &bdryNodeList, 
                   const std::vector<int> &innerNodeList, const int N ){

    // Preallocate Poisson-matrix 
    Eigen::SparseMatrix<dType> A(N*N, N*N);
    
    buildPoissonMatrix<mType, dType>(A, bdryNodeList, innerNodeList, N);
    
    // To do: how to deal with float-z,b? Then this map does not work, is there a generic way?
    // probably need to use Eigen::Matrix<T, N*N, 1> or so
//     Eigen::Matrix<T, 1, 25> zE;
//     Eigen::VectorXd zE(N*N);
//     Eigen::VectorXd bE = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
//     std::cout << bE << std::endl;
    
    // Now solve the Poisson equation using sparse Cholesky factorization
    // for later: solveWithGuess()
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<dType> > solver;
    zE = solver.compute(A).solve(bE);
    
    std::cout << "Solution: " << std::endl;
    std::cout << zE << std::endl;
    
    // Copy zE back to z
//     std::vector<T> z(N*N);
//     for(int i=0; i<N*N; i++)
//         z[i] = zE[i]; // Not too happy with this ^^
//    
//     return z;
}

// ################################################################################################


// get differentials by FD
// An idea might be to outsource the stencils to own functions, but I doubt it would 
// increase readability and/or performance
template <typename mType, typename dType> 
void minSurfOperator( Eigen::MatrixBase<mType> &outVec, const Eigen::MatrixBase<mType> &inVec, 
                      std::vector<int> innerNodeList, const int N ){ // I am really not married to the name of this function
//    std::vector<T> outVec(N*N);
   const dType h = 1./N;
   dType tmp = 0;
   
   
   // Something does not work in the following loop...
   
   for(auto& i: innerNodeList) {
       for(auto& j: innerNodeList) {
           std::cout << "\nheeelloooo\n";
           // tmp = (1+z_x^2)*z_yy
           tmp = (1 + pow((inVec[i+1+j*N] - inVec[i-1+j*N]) / (2*h), 2))
                 * (inVec[i+1+j*N] -2*inVec[i  +j*N]+ inVec[i-1+j*N]) / (h*h);
           // tmp -= 2*z_x*z_y*z_xy
           tmp -= 2* (inVec[i+1+ j   *N] - inVec[i-1+ j   *N]) / (2*h) // z_x
                   + (inVec[i  +(j+1)*N] - inVec[i  +(j-1)*N]) / (2*h) // z_y
                   *   (inVec[i+1+(j+1)*N] + inVec[i-1+(j-1)*N]\
                      - inVec[i+1+(j-1)*N] - inVec[i+1+(j+1)*N]) / (4*h); // z_xy
           // tmp += (1+z_y^2)*z_xx
           tmp += (1 + pow((inVec[i+(j+1)*N] - inVec[i+(j-1)*N]) / (2*h), 2))
                 * (inVec[i+(j+1)*N] -2*inVec[i+j*N]+ inVec[i+(j-1)*N]) / (h*h);
           std::cout << "tmp for " << i << ", " << j << " is " << tmp << std::endl;
           if (tmp > EPS)
               outVec[i+j*N] = tmp;
           else
               outVec[i+j*N] = 0;
       }
    }

    
//    return outVec;
}
    
    
// ################################################################################################
// implement (1+z_x^2)z_yy - ...
/*
template <typename T>
void runNewton( std::vector<T> &old, std::vector<T> &fresh ) {

    intermediate = minSurfOperator(old); // == residual
    if(intermediate.norm() < TOL) { 
        compute " tangent \ intermediate " := delta;
        fresh = old - delta; 
    }
    else
        ABORT;
    
}
*/
