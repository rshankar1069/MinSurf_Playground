#include<iostream>
#include<stdio.h>
#include<vector>
#include<math.h>
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


// Boundary nodes by lexicographical ordering
// To be extended for a) different boundaries
//                    b) MPI?
void setBdryNodes( std::vector<int> &bdryNodeList, const int N ) {
    int i=0, j=0;
    
    // First row
    for(i=0; i<N; i++) {
        bdryNodeList.push_back(i);
//         std::cout << "me: " << i << std::endl;
    }
    
    // outer edges
    j = 0;
    while (j<N*(N-1)) {
        j+=N;
        bdryNodeList.push_back(j);
        bdryNodeList.push_back(j+N-1);
//         std::cout << "me: " << j << std::endl;
    }
    // Last row
    for(i=0; i<N; i++) {
        bdryNodeList.push_back(j + i);
//         std::cout << "me: " << i+j << std::endl;
    }
}

// Inner nodes by lexicographical ordering
void setInnerNodes( std::vector<int> &innerNodeList, const int N ) {
    
    for(int i=1; i<N-1; i++) {
        for(int j=1; j<N-1; j++) {
            innerNodeList.push_back(i+j*N);
        }
    }
        
}

// Function to apply boundary conditions - needs to be extended -> Sankar
template <typename T>
void applyBC( const int BCtype, std::vector<T> &V , const int N, std::vector<int> &bdryNodeList ) {
    for(auto& i: bdryNodeList) {
        V.at(i) = sin( (double)i /N);
    }

}

template <typename T>
void buildPoissonMatrix( Eigen::SparseMatrix<T> &A, std::vector<int> &bdryNodeList, int N ) {
    typedef Eigen::Triplet<double> triplet;

	std::vector<triplet> tripletList;
	tripletList.reserve(3*N);
    double v_ij;

    // Assemble the FD matrix (should work, validated with matlab)
    // Set only inner nodes
    // NEED TO CHECK: NOW, DUE TO THE TRIPLETS, BOUNDARY ENTRIES ARE ADDED
    // ON TOP, WHICH SHOULD NOT BE THE CASE. WE ONLY WANT THEM TO BE 1
    // THIS NEEDS TO BE CORRECTED, MAYBE DO NOT USE TRIPLETS AFTER ALL
    // OR ACTUALLY USE THE innerNodeList!
	for(int i=1; i<(N-1)*(N-1); i++) {
        for(int j=1; j<(N-1)*(N-1); j++) {
            if (i==j)
                v_ij = 4;
            else if (i-1==j and i-1>=0)
                v_ij = -1;
            else if (i+1==j and i+1<N*N)
                v_ij = -1;
            else if (i+N-1==j and i+N-1<N*N)
                v_ij = -1;
            else if (i-N+1==j and i-N+1>=0)
                v_ij = -1;
            else v_ij = 0;
            
//             std::cout << i << "," << j << ": " << v_ij << std::endl;
            if (v_ij)
                tripletList.push_back(triplet(i,j,v_ij));
        }
    }
    // Set a 1 where Dirichlet BC applies
    for(auto& i: bdryNodeList)
        tripletList.push_back(triplet(i, i,1));

    // Build sparse A from triplets
    A.setFromTriplets(tripletList.begin(), tripletList.end());

}

template <typename T>
void getInitGuess( std::vector<T> &z, std::vector<T> &b, std::vector<int> &bdryNodeList, int N ){
    // Preallocate Poisson-matrix 
    Eigen::SparseMatrix<T> A(N*N, N*N);
    
    buildPoissonMatrix(A, bdryNodeList, N);
     
    Eigen::VectorXd zE = Eigen::Map<Eigen::VectorXd>(z.data(), z.size());
    Eigen::VectorXd bE = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
//     std::cout << bE << std::endl;
    
    // Now solve the Poisson equation using sparse Cholesky factorization
    // for later: solveWithGuess()
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    zE = solver.compute(A).solve(bE);
    
    std::cout << A << std::endl;
    
    // Copy zE back to z
    
    
}

/*
template <typename T>
T getDx...
// get differentials by FD
//

template <typename T> 
void minSurfOperator(){}
// implement (1+z_x^2)z_yy - ...

*/
