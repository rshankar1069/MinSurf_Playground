/*
 *
 *
 *
 */

#include "minSurf.h"


// Whatever only needs to touch the inner nodes, maybe (and probably)
//       there is a smart way to avoid allocating the full N*N
// Returning the vector of the minSurfOperator does not really work
// --> Switch to Eigen-containers!
//     --> Now computing the minSurfOperator explodes, but should be resolved quickly


int main() {
    
    // Try to solve Poisson equation using Eigen
    const int N = 5; // 100 within 1sec for Poisson, but 500 intractable...
    typedef double dType;
    typedef Eigen::Matrix<dType, N*N, 1> Vector;
    std::vector<int> bdryNodeList, innerNodeList;
    setBdryNodes(bdryNodeList, N);
    setInnerNodes(innerNodeList, N);
    
    
    
    // Prepare solution vector and RHS with BC
//     std::vector<double> z(N*N), b(N*N);
    Vector z, b;
    applyBC(0, b, innerNodeList, bdryNodeList);
    
    getInitGuess<Vector, dType>(z, b, bdryNodeList, innerNodeList, N);
   
    std::cout << "Initial guess is: " << std::endl;
//     for(auto& elem: z)
//         std::cout << "\t" << elem << std::endl;
    std::cout << z << std::endl;
    
    // Apply discrete minSurf on initial guess
    Vector rh;
    minSurfOperator<Vector, dType>(rh, z, innerNodeList, N);
    std::cout << "Differential is: " << std::endl;
//     for(auto& elem: rh)
//         std::cout << "\t" << elem << std::endl;
    std::cout << rh << std::endl;
   
    // Run Newtons method
    //solve(...);
}
