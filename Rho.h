#ifndef RHO_H
#define RHO_H

#include "Basis.h"
#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Naive density reconstruction helper.
 *
 * Provides a straight implementation of the six nested loops described in the
 * project brief to accumulate the local density from the density matrix.
 */
class Rho
{
public:
    /**
     * @brief Compute the density on (r,z) grids using the naive algorithm.
     * @param basis Truncated basis definition with quantum number limits.
     * @param rho Density matrix of size (basis size) x (basis size).
     * @param rVals Sampling points along the radial axis.
     * @param zVals Sampling points along the axial axis.
     * @return Matrix of shape (rVals.n_elem, zVals.n_elem) containing ρ(r,z).
     */
    static arma::mat compute(const Basis& basis,
                             const arma::mat& rho,
                             const arma::vec& rVals,
                             const arma::vec& zVals);
};

#endif // RHO_H
