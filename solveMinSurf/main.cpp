/*  FYI, for now, I just compiled with g++ -std=c++14 -Wall
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

    // Set boundary and inner nodes (based on structure grid 
    // and lexicographical ordering
    typedef std::vector<int> listType;
    listType bdryNodeList, innerNodeList;
    setBdryNodes(bdryNodeList, N);
    setInnerNodes(innerNodeList, N);
    
    
    
    // Prepare solution vector and RHS with BC
    Vector z=Vector::Zero();
    Vector b = Vector::Zero();
    applyBC<Vector, dType, listType>(0, b, innerNodeList, bdryNodeList);
    
    getInitGuess<Vector, dType, listType>(z, b, bdryNodeList, innerNodeList, N);
   
//     std::cout << "Initial guess is: " << std::endl;
//     std::cout << z << std::endl;
    
    // Apply discrete minSurf on initial guess
    Vector rh = Vector::Zero();
    minSurfOperator<Vector, dType, listType>(rh, z, innerNodeList, N);
//     std::cout << "Differential is: " << std::endl;
//     std::cout << rh << std::endl;
    
    
    // Run Newtons method
    //solve(...);
}
