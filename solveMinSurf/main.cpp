/*
 *
 *
 *
 */

#include "minSurf.h"


int main() {
	
	// Try to solve Poisson equation using Eigen
	int N = 5; // 100 within 1sec for Poisson, but 500 intractable...
    std::vector<int> bdryNodeList, innerNodeList;
	setBdryNodes(bdryNodeList, N);
    setInnerNodes(innerNodeList, N);
    
    
    
    // Prepare solution vector and RHS with BC
    std::vector<double> z(N*N), b(N*N);
    std::fill(z.begin(), z.end(), 0);
    std::fill(b.begin(), b.end(), 0);
    applyBC(0, b, N, bdryNodeList);
    
    getInitGuess(z, b, bdryNodeList, N);
    
    
   
    std::cout << "Solution is: " << std::endl;
    for(auto& elem: z)
        std::cout << "\t" << elem << std::endl;
    
}
