#include<iostream>
#include<stdio.h>
#include<vector>
#include<valarray>
#include<cmath>
#include<fstream>
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/Core"
#include"Eigen/Eigen/SparseCholesky"


#define EPS 1.e-12


template <typename T>
T residual( const std::vector<T> &z, std::vector<T> &r ) {
    // computes residual entries in r[i]
    // returns sum of squares of r
    T sum = 0;
    for (auto& n : r)
        sum += n*n;
}

// ################################################################################################
// Boundary nodes by lexicographical ordering
// To be extended for a) different boundaries
//                    b) MPI?
template <typename listType>
void setBdryNodes( listType &bdryNodeList, const int N ) {
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
    
    // Sort (works for vector) -> why? my hope is that data is better aligned in the following then
    // and it is worth the effort
    std::sort(bdryNodeList.begin(), bdryNodeList.end());

    // Debug output
    for(auto& i: bdryNodeList)
        std::cout << i << std::endl;
}

// Inner nodes by lexicographical ordering
template <typename listType>
void setInnerNodes( listType &innerNodeList, const int N ) {
    
    // Set inner nodes (assuming structured cartesian grid)
    for(int i=1; i<N-1; i++) {
        for(int j=1; j<N-1; j++) {
            innerNodeList.push_back(i+j*N);
        }
    }
    // Sort list (works for vector)
    std::sort(innerNodeList.begin(), innerNodeList.end());
    // Debug output
    std::cout << "innernodes:\n";
    for(auto& it: innerNodeList)
        std::cout << it << std::endl;
        
}

// Function to apply boundary conditions - needs to be extended -> Sankar
template <typename mType, typename dType, typename listType>
void applyBC( const int BCtype, Eigen::MatrixBase<mType> &V , const listType &bdryNodeList, 
            listType innerNodeList ) {
    // Set boundary values
    for(auto& i: bdryNodeList)
        V[i]= 5;//sin( (dType)i );
    // Set rest to zero
    for(auto& i: innerNodeList)
        V[i] = 0;
}

// ################################################################################################
template <typename mType, typename dType, typename listType>
void buildPoissonMatrix( Eigen::SparseMatrix<dType> &A, const listType &bdryNodeList, 
                         const listType &innerNodeList, const int N ) { 

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

// Set a 1 where Dirichlet BC applies
    for(auto& i: bdryNodeList)
        tripletList.push_back(triplet(i, i,1));

    // Build sparse A from triplets
    A.setFromTriplets(tripletList.begin(), tripletList.end());
//     std::cout << A << std::endl;
}

template <typename mType, typename dType, typename listType>
void getInitGuess( Eigen::MatrixBase<mType> &zE, const Eigen::MatrixBase<mType> &bE, const listType &bdryNodeList, 
                   const listType &innerNodeList, const int N ){

    // Preallocate Poisson-matrix 
    Eigen::SparseMatrix<dType> A(N*N, N*N);
    
    buildPoissonMatrix<mType, dType, listType>(A, bdryNodeList, innerNodeList, N);
        
    // Now solve the Poisson equation using sparse Cholesky factorization
    // for later: solveWithGuess()
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<dType> > solver;
    zE = solver.compute(A).solve(bE);
    
}

// ################################################################################################


// Get differentials by FD
// An idea might be to outsource the stencils to own functions, but I doubt it would 
// increase readability and/or performance
template <typename mType, typename dType, typename listType> 
void minSurfOperator( Eigen::MatrixBase<mType> &outVec, const Eigen::MatrixBase<mType> &inVec, 
                      const listType innerNodeList, const int N ){ // I am really not married to the name of this function
   const dType h = 1./N;
   dType tmp = 0;
   
   // Maybe we can exploit Eigen a little bit more to make the index-accessing a little more
   // convenient or more efficient...
   
   // Also, the result (if correct, but I think so...) is very sparse!
   // Does it maybe make sense to store it in terms of "duplets" -> tuple of (index, value)??

   for(auto& i: innerNodeList) {
       // tmp = (1+z_x^2)*z_yy
       tmp = (1 + pow((inVec[i+1] - inVec[i-1]) / (2*h), 2))
             * (inVec[i+N] -2*inVec[i]+ inVec[i-N]) / (h*h);
       // tmp -= 2*z_x*z_y*z_xy
       tmp -= 2* (inVec[i+1] - inVec[i-1]) / (2*h) // z_x
               * (inVec[i+N] - inVec[i-N]) / (2*h) // z_y
               *   (inVec[i+1+N] + inVec[i-1-N]
                  - inVec[i+1-N] - inVec[i-1+N]) / (4*h); // z_xy
       // tmp += (1+z_y^2)*z_xx
       tmp += (1 + pow((inVec[i+N] - inVec[i-N]) / (2*h), 2))
            * (inVec[i+1] -2*inVec[i]+ inVec[i-1]) / (h*h);
       if (fabs(tmp) > EPS)
           outVec[i] = tmp;
       else
           outVec[i] = 0.;
    }
    
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
