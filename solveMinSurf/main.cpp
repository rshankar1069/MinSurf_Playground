/*
 *
 *
 *
 */

#include "minSurf.h"


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
    std::valarray<double> rh = minSurfOperator(z, innerNodeList, N);
    for(auto& elem: rh)
        std::cout << "\t" << elem << std::endl;
    
   
    // Run Newtons method
    //solve(...);
}
