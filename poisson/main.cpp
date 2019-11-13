/*
 *
 *
 *
 */

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


// Maybe make a template 
// template <class T>
// class solutionVector


// Boundary nodes by lexicographical ordering
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

// Function to apply boundary conditions - needs to be extended
template <typename T>
void applyBC( const int BCtype, std::vector<T> &V , const int N, std::vector<int> &bdryNodeList ) {
    for(int i=0; i<bdryNodeList.size(); i++) {
        V.at(bdryNodeList.at(i)) = sin( (double)bdryNodeList.at(i) /N);
    }

}


int main() {
	
	// Try to solve Poisson equation using Eigen
	int N = 100; // 100 within 1sec for Poisson, but 500 intractable...
    std::vector<int> bdryNodeList;
	setBdryNodes(bdryNodeList, N);
    typedef Eigen::Triplet<double> triplet;

	std::vector<triplet> tripletList;
	tripletList.reserve(3*N);
    double v_ij;

    // Assemble the FD matrix (should work, validated with matlab)
    // Might be done more efficiently though...
	for(int i=0; i<N*N; i++) {
        for(int j=0; j<N*N; j++) {
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

    
    Eigen::SparseMatrix<double> A(N*N, N*N);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
//     std::cout << "Matrix: " << A << std::endl;
    
    // Prepare solution vector and RHS with BC
    std::vector<double> z(N*N), b(N*N);
    std::fill(z.begin(), z.end(), 0);
    std::fill(b.begin(), b.end(), 0);
    applyBC(0, b, N, bdryNodeList);
    
    Eigen::VectorXd zE = Eigen::Map<Eigen::VectorXd>(z.data(), z.size());
    Eigen::VectorXd bE = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
//     std::cout << bE << std::endl;
    
    // Now solve the Poisson equation using sparse Cholesky factorization
    // for later: solveWithGuess()
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    zE = solver.compute(A).solve(bE);
   
    std::cout << "Solution is: " << std::endl << zE << std::endl;
}
