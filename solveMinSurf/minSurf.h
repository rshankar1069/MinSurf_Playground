#include<iostream>
#include<stdio.h>
#include<vector>
#include<valarray>
#include<cmath>
#include<fstream>
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/Core"
#include"Eigen/Eigen/SparseCholesky"


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
template <typename T>
void applyBC( const int BCtype, std::vector<T> &V , std::vector<int> bdryNodeList, 
            std::vector<int> innerNodeList ) {
    // Set boundary values
    for(auto& i: bdryNodeList)
        V[i]= sin( (double)i );
    // Set rest to zero
    for(auto& i: innerNodeList)
        V[i] = 0;
}

// ################################################################################################
template <typename T>
void buildPoissonMatrix( Eigen::SparseMatrix<T> &A, std::vector<int> bdryNodeList, 
                         std::vector<int> innerNodeList, const int N ) {

    typedef Eigen::Triplet<double> triplet;
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

template <typename T>
std::vector<T> getInitGuess( std::vector<T> b, std::vector<int> bdryNodeList, 
                   std::vector<int> innerNodeList, const int N ){

    // Preallocate Poisson-matrix 
    Eigen::SparseMatrix<T> A(N*N, N*N);
    
    buildPoissonMatrix(A, bdryNodeList, innerNodeList, N);
    
    // To do: how to deal with float-z,b? Then this map does not work, is there a generic way?
    // probably need to use Eigen::Matrix<T, N*N, 1> or so
    Eigen::Matrix<T, 1, 25> zE;
//     Eigen::VectorXd zE(N*N);
    Eigen::VectorXd bE = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
//     std::cout << bE << std::endl;
    
    // Now solve the Poisson equation using sparse Cholesky factorization
    // for later: solveWithGuess()
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    zE = solver.compute(A).solve(bE);
    
    std::cout << "Solution: " << std::endl;
    std::cout << zE << std::endl;
    
    // Copy zE back to z
    std::vector<T> z(N*N);
    for(int i=0; i<N*N; i++)
        z[i] = zE[i]; // Not too happy with this ^^
//    
    return z;
}

// ################################################################################################
template <typename T>
std::vector<T> getDxP( std::vector<T> z, std::vector<int> innerNodeList, const double h, 
                       const int N ) {
    std::vector<T> Dx(N*N);
    
    for(auto& i: innerNodeList) {
        Dx[i] = (z[i+1] - z[i]) / (2*h);
    }
    
    return Dx;
}

template <typename T>
std::vector<T> getDxP( std::vector<T> z, std::vector<int> innerNodeList, const double h, 
                       const int N, const int power ) {
    std::vector<T> Dx(N*N);
    
    for(auto& i: innerNodeList) {
        Dx[i] = pow((z[i+1] - z[i]) / (2*h), power);
    }
    
    return Dx;
}

// Without power, i.e z_y
template <typename T>
std::valarray<T> getDyP( std::vector<T> z, std::vector<int> innerNodeList, const double h, 
                         const int N ) {
    std::valarray<T> Dy(N*N);
    
    for(auto& i: innerNodeList) {
        Dy[i] = (z[i+N] - z[i-N]) / (2*h);
    }
    
    return Dy;
}

// With power, i.e. useful for z_y^2
template <typename T>
std::valarray<T> getDyP( std::vector<T> z, std::vector<int> innerNodeList, const double h, 
                         const int N, const int power ) {
    std::valarray<T> Dy(N*N);
    
    for(auto& i: innerNodeList) {
        Dy[i] = pow((z[i+N] - z[i-N]) / (2*h), power);
    }
    
    return Dy;
}

template <typename T>
std::valarray<T> getDxx( std::vector<T> z, std::vector<int> innerNodeList, const double h, 
                         const int N ) {
    std::valarray<T> Dxx(N*N);
    
    for(auto& i: innerNodeList) {
        Dxx[i] = (z[i+1] - 2*z[i] + z[i]) / (h*h);
    }
    
    return Dxx;
    
}

template <typename T>
std::valarray<T> getDyy( std::vector<T> z, std::vector<int> innerNodeList, const double h, 
                         const int N ) {
    std::valarray<T> Dyy(N*N);
    
    for(auto& i: innerNodeList) {
        Dyy[i] = (z[i+N] - 2*z[i] + z[i-N]) / (h*h);
    }
    
    return Dyy;
    
}


// get differentials by FD
//
// Seriously need to rewrite this...
template <typename T> 
std::valarray<T> minSurfOperator( std::vector<T> inVec, std::vector<int> innerNodeList, const int N ){ // I am really not married to the name of this function
   std::valarray<T> outVec;
   const double h = 1./N;
   
   outVec = (1.+getDxP(inVec, innerNodeList, h, N, 2))*getDyy(inVec, innerNodeList, h, N) 
          - 2.*getDxP(inVec, innerNodeList, h, N)*getDy(inVec, innerNodeList, h, N)*getDxDy(inVec, innerNodeList, h, N)
          + (1+getDy2(inVec, innerNodeList, h, N, 2))*getDxx(inVec, innerNodeList, h, N);
          
    // If you don't want to use the valarrays to continue here, maybe convert/copy everything now
    // But I think it makes sense to use them, since we can directly add/multiply vectors here.
    // 2 points to check: -> Does the code above really do what we want (not sure about the (1+z_x^2)*z_yy in the vector mult. notation
    //                    -> Do we maybe use Eigen-datatypes in general? these could add/substract/multiply aswell...
 
    return outVec;
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
