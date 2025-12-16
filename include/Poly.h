#ifndef POLY_H
#define POLY_H

#include "../armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Utility class for generating Hermite and Laguerre polynomial tables.
 *
 * The class precomputes polynomial values on supplied grids so subsequent
 * lookups are O(1) without redundant recurrence evaluations.
 */
class Poly
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::mat hermiteTable; 
    arma::cube laguerreTable;

    /**
     * @brief Compute Hermite polynomials up to order n on the provided points.
     * @param n highest order to tabulate (inclusive).
     * @param zVals vector of evaluation points.
     */
    void calcHermite(int n , const arma::vec& zVals);

    /**
     * @brief Compute generalized Laguerre polynomials for multiple (m,n) pairs.
     * @param mCount number of alpha values (rows) to generate.
     * @param nCount number of polynomial orders per alpha.
     * @param zVals vector of evaluation points.
     */
    void calcLaguerre(int mCount , int nCount , const arma::vec& zVals);

    /**
     * @brief Debug helper that prints a dense matrix to stdout.
     */
    void static printMatrix(arma::mat mat);

    /**
     * @brief Placeholder recurrence hook for custom Laguerre slice building.
     */
    arma::mat calcSliceN(int n , arma::mat sliceNMoins1 , arma::mat sliceNMoins2); 

    /**
     * @brief Retrieve the cached Hermite polynomial of the requested order.
     */
    arma::vec hermite(int order) const;

    /**
     * @brief Return the cached Hermite row as a column vector.
     */
    arma::vec getHermiteRow(int order) const;

    /**
     * @brief Retrieve the cached Laguerre tube for the (m,n) pair.
     */
    arma::vec laguerre(int mOrder, int nOrder) const;


};

#endif
