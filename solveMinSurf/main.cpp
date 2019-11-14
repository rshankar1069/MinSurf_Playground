/*
 *
 *
 *
 */

#include "minSurf.h"


// To do: Use Eigen-containers for the numerical important values, i.e. b,z,r,...
// Whatever only needs to touch the inner nodes, maybe (and probably)
//       there is a smart way to avoid allocating the full N*N
// Returning the vector of the minSurfOperator does not really work
// --> Switch to Eigen-containers!

int main() {
    
    // Try to solve Poisson equation using Eigen
    const int N = 5; // 100 within 1sec for Poisson, but 500 intractable...
    std::vector<int> bdryNodeList, innerNodeList;
    setBdryNodes(bdryNodeList, N);
    setInnerNodes(innerNodeList, N);
    
    
    
    // Prepare solution vector and RHS with BC
    std::vector<double> z(N*N), b(N*N);
    applyBC(0, b, innerNodeList, bdryNodeList);
    
    z = getInitGuess(b, bdryNodeList, innerNodeList, N);
   
    std::cout << "Initial guess is: " << std::endl;
    for(auto& elem: z)
        std::cout << "\t" << elem << std::endl;
    
    // Apply discrete minSurf on initial guess
    std::vector<double> rh = minSurfOperator(z, innerNodeList, N);
    std::cout << "Differential is: " << std::endl;
    for(auto& elem: rh)
        std::cout << "\t" << elem << std::endl;
    
   
    // Run Newtons method
    //solve(...);
}
