#include <armadillo>
#include "solver.h"
#include "basis.h"



/**
 * @brief Class constructor.
 * It loads the rho matrix, and fills the indices matrix.
 * @param _basis Instance of the basis to use.
 */
Solver::Solver(Basis _basis): basis(_basis) {
    // load the rho file and place it inside rho matrix
    rho.load("input/rho.arma", arma::arma_ascii);
    indices = arma::icube(basis.mMax, basis.nMax.max(), basis.n_zMax.max());
    uint i = 0;
    for (int m = 0; m < basis.mMax; m++)
        for (int n = 0; n < basis.nMax(m); n++)
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                indices(m, n, n_z) = i;
                i++;
            }
}

/**
 * @brief Resets the variables tracking the maximum degrees of polynomials that has been calculated.
 * This function must be called if you wan't to rerun a version function with different values for r and z so the call to the version function won´t try to use already calculated values for Laguerre and hermite polynomials .
 */
void Solver::reset() {
    basis.calc_up_to_m  = 0;
    basis.calc_up_to_n  = 0;
    basis.calc_up_to_nz = 0;
}


/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This is the first version, a naive algorithm.
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version0(arma::Col<double> rVals, arma::Col<double> zVals) {
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    // a and b are used to convert
    uint a = 0;
    uint b = 0;
    // loop over all indexes
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                for (int mp = 0; mp < basis.mMax; mp++) {
                    for (int np = 0; np < basis.nMax(mp); np++) {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++) {
                            // calculate density for those indexes
                            arma::mat funcA = basis.basisFunc( m,  n,  n_z, zVals, rVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            result += funcA % funcB * rho(a, b); // mat += mat % mat * double
                            b++;
                        } // n_zp
                    } // np
                } // mp
                b=0;
                a++;
            } // n_z
        } // n
    } // m
    return result;
}
