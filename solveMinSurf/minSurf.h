#include <iostream>
#include <stdio.h>
#include <vector>
#include <valarray>
#include <cmath>
#include <fstream>
#include "Eigen/Eigen/Sparse"
#include "Eigen/Eigen/Eigenvalues"
#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/SparseLU"
#include "Eigen/Eigen/SparseCholesky"
#include "Eigen/Eigen/IterativeLinearSolvers"

#define EPS 1.e-12

// #################################################################################################
// Boundary nodes by lexicographical ordering
// To be extended for a) different boundaries
//                    b) MPI?
template <typename listType>
void setBdryNodes(listType &bdryNodeList, const int N)
{
    int i = 0;

    // Set first and last row
    for (i = 0; i < N; i++)
    {
        // First row
        bdryNodeList.push_back(i);
        // Last row
        bdryNodeList.push_back(N * (N - 1) + i);
    }

    // outer edges
    for (i = N; i < N * (N - 1); i += N)
    {
        bdryNodeList.push_back(i);
        bdryNodeList.push_back(i + N - 1);
    }

    // Sort (works for vector)
    std::sort(bdryNodeList.begin(), bdryNodeList.end());
}

// Inner nodes by lexicographical ordering
template <typename listType>
void setInnerNodes(listType &innerNodeList, const int N)
{

    // Set inner nodes (assuming structured cartesian grid)
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < N - 1; j++)
        {
            innerNodeList.push_back(i + j * N);
        }
    }
    // Sort list
    std::sort(innerNodeList.begin(), innerNodeList.end());
}

// Function to apply boundary conditions
template <typename mType, typename dType, typename listType>
void applyBC(const int BCtype, Eigen::MatrixBase<mType> &V, const listType &innerNodeList,
             listType bdryNodeList, const int N)
{
    // Set boundary values
    int i,j; dType x, y;
    for(auto& node: bdryNodeList) {
        i = node%N;
        x = ((dType) i) / (N-1);
        j = node/N;
        y = ((dType) j) / (N-1);
        
        if( i==0 or i==(N-1)) // lower or upper bdry -> x^2
            V[node] = -.1*pow(y-.5,2);
        else if ( j==0 or j==(N-1)) // left or right bdry -> y^2
            V[node] = -.1*pow(x-.5,2);
    }
}

// ################################################################################################
template <typename mType, typename dType, typename listType>
void buildPoissonMatrix(Eigen::SparseMatrix<dType> &A, const listType &innerNodeList,
                        const listType &bdryNodeList, const int N)
{

    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(5 * N);

    // Assemble the FD matrix
    // Set inner nodes - 5-pt stencil of FD
    for (auto &i : innerNodeList)
    {
        tripletList.push_back(triplet(i, i, 4));
        tripletList.push_back(triplet(i, i + 1, -1));
        tripletList.push_back(triplet(i, i - 1, -1));
        tripletList.push_back(triplet(i, i + N, -1));
        tripletList.push_back(triplet(i, i - N, -1));
    }

    // Set a 1 where Dirichlet BC applies
    for (auto &i : bdryNodeList)
        tripletList.push_back(triplet(i, i, 1));

    // Build sparse A from triplets
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <typename mType, typename dType, typename listType>
void getInitGuess(Eigen::MatrixBase<mType> &zE, const Eigen::MatrixBase<mType> &bE,
                  const listType &innerNodeList, const listType &bdryNodeList, const int N)
{

    // Preallocate Poisson-matrix
    Eigen::SparseMatrix<dType> A(N * N, N * N);

    buildPoissonMatrix<mType, dType, listType>(A, innerNodeList, bdryNodeList, N);

    // Now solve the Poisson equation using sparse Cholesky factorization
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<dType>> solver;
    zE = solver.compute(A).solve(bE);
}

// ################################################################################################

// Action of discrete minSurfOperator on an input vector
template <typename mType, typename dType, typename listType>
void minSurfOperator(Eigen::MatrixBase<mType> &y, const Eigen::MatrixBase<mType> &x, Eigen::SparseMatrix<dType> &Jacobian,
                     const listType innerNodeList, const int N)
{
    const dType h = 1. / ((dType)N);

    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(9 * N * N);

    dType a1_y;
    std::vector<dType> a1_x(9);

    for (auto &i : innerNodeList)
    {
        // forward
        // d(...) = fd(x)
        dType dx = (x[i + 1] - x[i - 1]) / (2 * h);
        dType dy = (x[i + N] - x[i - N]) / (2 * h);
        dType dxx = (x[i + 1] - 2 * x[i] + x[i - 1]) / (h * h);
        dType dyy = (x[i + N] - 2 * x[i] + x[i - N]) / (h * h);
        dType dxy = (x[i + 1 + N] + x[i - 1 - N] - x[i + 1 - N] - x[i - 1 + N]) / (4 * h * h);

        // v1 = (1+z_x^2)*z_yy
        dType v1 = (1 + pow(dx, 2)) * dyy;
        // v2 = -2*z_x*z_y*z_xy
        dType v2 = -2 * dx * dy * dxy;
        // v3 = (1+z_y^2)*z_xx
        dType v3 = (1 + pow(dy, 2)) * dxx;

        y[i] = v1 + v2 + v3;

        // reverse
        // seed a1_y
        a1_y = 1.;
    
        // reverse of y[i] = v1 + v2 + v3
        dType a1_v1 = a1_y;
        dType a1_v2 = a1_y;
        dType a1_v3 = a1_y;

        //reverse of v3 = (1 + pow(dy, 2)) * dxx
        dType a1_dy = dxx * 2 * dy * a1_v3;
        dType a1_dxx = (1 + pow(dy, 2)) * a1_v3;

        //reverse v2 = -2 * dx * dy * dxy
        dType a1_dx = -2 * dy * dxy * a1_v2;
        a1_dy += -2 * dx * dxy * a1_v2;
        dType a1_dxy = -2 * dx * dy * a1_v2;

        //reverse of v1 = (1 + pow(dx, 2)) * dyy
        a1_dx += dyy * 2 * dx * a1_v1;
        dType a1_dyy = (1 + pow(dx, 2)) * a1_v1;

        // reverse of d(...) = fd(x)
        a1_x[5] = a1_dx / (2 * h);
        a1_x[3] = a1_dx / (-2 * h);
        a1_x[7] = a1_dy / (2 * h);
        a1_x[1] = a1_dy / (-2 * h);
        a1_x[5] += a1_dxx / (h * h);
        a1_x[3] += a1_dxx / (h * h);
        a1_x[4] = (-2) * a1_dxx / (h * h);
        a1_x[7] += a1_dyy / (h * h);
        a1_x[1] += a1_dyy / (h * h);
        a1_x[4] += (-2) * a1_dyy / (h * h);
        a1_x[8] = a1_dxy / (4 * h * h);
        a1_x[0] = a1_dxy / (4 * h * h);
        a1_x[2] = a1_dxy / (-4 * h * h);
        a1_x[6] = a1_dxy / (-4 * h * h);

        // store derivatives in triplets
        tripletList.push_back(triplet(i, i - 1 - N, a1_x[0]));
        tripletList.push_back(triplet(i, i     - N, a1_x[1]));
        tripletList.push_back(triplet(i, i + 1 - N, a1_x[2]));
        tripletList.push_back(triplet(i, i - 1,     a1_x[3]));
        tripletList.push_back(triplet(i, i,         a1_x[4]));
        tripletList.push_back(triplet(i, i + 1,     a1_x[5]));
        tripletList.push_back(triplet(i, i - 1 + N, a1_x[6]));
        tripletList.push_back(triplet(i, i     + N, a1_x[7]));
        tripletList.push_back(triplet(i, i + 1 + N, a1_x[8]));
    }
    // Build sparse matrix from triplets
    Jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
}

// #################################################################################################
// minSurf solving routine

// Get residual and Jacobian by application of minSurf-operator
template <typename mType, typename dType, typename listType>
dType residual(Eigen::MatrixBase<mType> &resVec, const Eigen::MatrixBase<mType> &solVec, Eigen::SparseMatrix<dType> &Jacobian, 
               const listType &innerNodeList, const int N)
{
    // computes residual entries in resVec and Jacobian
    // returns norm of r

    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and
    // the boundary information on the bdryNodeList
    minSurfOperator<mType, dType, listType>(resVec, solVec, Jacobian, innerNodeList, N);

    // maybe one could out-source this applyBC, since it is not touched again...
    // but then, no setting resVec to zero and careful with setZero in minSurfOperator
    dType res = 0;
    for (auto &i : innerNodeList)
        res += pow(resVec[i], 2);

    return sqrt(res);
}

// Main solver loop
template <typename mType, typename dType, typename listType>
void runSolver(const listType &innerNodeList, const listType &bdryNodeList, const int N)
{
    mType z = mType::Zero();
    mType b = mType::Zero();

    applyBC<mType, dType, listType>(0, b, innerNodeList, bdryNodeList, N);

    getInitGuess<mType, dType, listType>(z, b, innerNodeList, bdryNodeList, N);

    applyBC<mType, dType, listType>(0, z, innerNodeList, bdryNodeList, N);

    dType res;
    mType resVec, dz;
    Eigen::SparseMatrix<dType> Jacobian(N * N, N * N);

    res = residual<mType, dType, listType>(resVec, z, Jacobian, innerNodeList, N);

    std::cout << "Starting residual: " << res << std::endl;

    dType omega = 0.85; // Relaxation parameter for Newton-Raphson
    unsigned iterationIndex = 0;
    // In case initiall guess was not horrifically lucky, run Newton-Raphson
    do
    {
        // Test for Eigenvalues of Jacobian - only test purpose, to know whether CG is a good idea or not

        // dz_n = grad[F(z_n)]^-1 * F(z_n)
        // To be played with: preconditioner (MUST),
        // initial guess (maybe inversion of the Poisson-gradient might also help, but no idea),
        //     tolerance (MUST).. should not be too high, as our main goal is the result of Newton
        Eigen::BiCGSTAB<Eigen::SparseMatrix<dType>> solver;
        solver.compute(Jacobian);
        dz = solver.solve(resVec);
        z -= omega * dz;
        applyBC<mType, dType, listType>(0, z, innerNodeList, bdryNodeList, N);

        // get residual, resVec and Jacobian -> F(z_n), grad[F(z_n)]
        res = residual<mType, dType, listType>(resVec, z, Jacobian, innerNodeList, N);

        iterationIndex++;
        if (!(iterationIndex % 100))
            std::cout << "\tAt iteration " << iterationIndex << " res is " << res << std::endl;
    } while (res > 1.e-6 && iterationIndex < 100);

    std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
              << res << "." << std::endl;
    std::cout << z << std::endl;
}
