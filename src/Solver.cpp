#include "../include/Solver.h"

/**
 * @brief Class constructor.
 * It loads the rho matrix, and fills the indices matrix.
 * @param _Basis Instance of the Basis to use.
 */
Solver::Solver(const Basis & _Basis) : basis(_Basis) {
    // load the rho file and place it inside rho matrix
    rho.load("rho.arma", arma::arma_ascii);
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
void Solver::reset() {}


/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This is the first version, a naive algorithm.
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version0(const arma::vec& rVals, const arma::vec& zVals) {
    arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem); // number of points on r- and z- axes
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



/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This version improves on the naive algorithm by leveraging the fact that \a rho is m-diagonal.
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version1(const arma::vec& rVals, const arma::vec& zVals) {
    arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem); // number of points on r- and z- axes
    // a and b are used to convert
    uint a = 0;
    uint b = 0;
    uint n_zMax_sum = 0;
    // loop over all indexes
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                b = n_zMax_sum;
                for (int np = 0; np < basis.nMax(m); np++) {
                    for (int n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++) {
                        // calculate density for those indexes
                        arma::mat funcA = basis.basisFunc(m, n, n_z, zVals, rVals);
                        arma::mat funcB = basis.basisFunc(m, np, n_zp, zVals, rVals); // mp is replaced by m because rho is m-diagonal
                        result += funcA % funcB * rho(a, b); // mat += mat % mat * double
                        b++;
                    } // n_zp
                } // np
                a++;
            } // n_z
        } // n
        n_zMax_sum += arma::sum(basis.n_zMax.row(m));
    } // m
    return result;
}



/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This version improves on the naive algorithm by using the fact that \a rho is both symmetric and m-diagonal.
 * It also uses a precalculated matrix of indexes to convert the unique set ( \a m, \a n, \a n_z) into a unique int,
 * instead of calculating it on the fly.
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version2(const arma::vec& rVals, const arma::vec& zVals) {
    arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem); // number of points on r- and z- axes
    // a and b are used to convert
    uint a = 0;
    uint b = 0;
    // loop over all indexes
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                for (int np = 0; np < n+1; np++) {
                    for (int n_zp = 0; n_zp < (np < n ? basis.n_zMax(m, np) : n_z+1); n_zp++) {
                        a = indices(m, n, n_z);
                        b = indices(m, np, n_zp);
                        // calculate density for those indexes
                        arma::mat funcA = basis.basisFunc( m,  n,  n_z, zVals, rVals);
                        arma::mat funcB = basis.basisFunc(m, np, n_zp, zVals, rVals);
                        if (a == b) {
                            result += funcA % funcB * rho(a, b); // mat += mat % mat * double
                        } else if (a > b) {
                            result += 2 * funcA % funcB * rho(a, b);
                        }
                    } // n_zp
                } // np
            } // n_z
        } // n
    } // m
    return result;
}



/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This version improve on version2 by calculating some values higher in the nested loops (and thus sooner).
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version3(const arma::vec& rVals, const arma::vec& zVals) {
    arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem); // number of points on r- and z- axes
    // a and b are used to convert
    uint a = 0;
    uint b = 0;
    // loop over all indexes
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                a = indices(m, n, n_z);
                arma::mat funcA = basis.basisFunc( m,  n,  n_z, zVals, rVals);
                for (int np = 0; np < n+1; np++) {
                    for (int n_zp = 0; n_zp < (np < n ? basis.n_zMax(m, np) : n_z+1); n_zp++) {
                        b = indices(m, np, n_zp);
                        // calculate density for those indexes
                        arma::mat funcB = basis.basisFunc(m, np, n_zp, zVals, rVals);
                        if (a == b) {
                            result += funcA % funcB * rho(a, b); // mat += mat % mat * double
                        } else if (a > b) {
                            result += 2 * funcA % funcB * rho(a, b);
                        }
                    } // n_zp
                } // np
            } // n_z
        } // n
    } // m
    return result;
}


/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This version improve on version3 by removing all of the conditional tests except the ones in function calls.
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version4(const arma::vec& rVals, const arma::vec& zVals) {
    arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem); // number of points on r- and z- axes
    // a and b are used to convert
    uint a = 0;
    uint b = 0;

    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                a = indices(m, n, n_z);
                arma::mat funcA = basis.basisFunc(m, n, n_z, zVals, rVals);
                result += funcA % funcA * rho(a, a); // diagonal
                for (int np = 0; np < n; np++) {
                    for (int n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++) {
                        b = indices(m, np, n_zp);
                        arma::mat funcB = basis.basisFunc(m, np, n_zp, zVals, rVals);
                        result += 2 * funcA % funcB * rho(a, b);
                    }
                }
                for (int n_zp = 0; n_zp < n_z; n_zp++) {
                    b = indices(m, n, n_zp);
                    arma::mat funcB = basis.basisFunc(m, n, n_zp, zVals, rVals);
                    result += 2 * funcA % funcB * rho(a, b);
                }
            }
        }
    }
    return result;
}


/**
 * @brief Calculates the density for each coordinates described by \a rVals and \a zVals.
 *
 * This version improve on version4 by not calculating moving calls to rPart higher up in the nested loops.
 *
 * @param rVals Vector of coordinates on the \a r axis
 * @param zVals Vector of coordinates on the \a z axis
 */
arma::mat Solver::version5(const arma::vec& rVals, const arma::vec& zVals) {
    arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem); // number of points on r- and z- axes
    // a and b are used to convert
    uint a = 0;
    uint b = 0;

    // variables to store calculated values
    arma::vec rPartA;
    arma::vec zPartA;
    arma::vec rPartB;
    arma::vec zPartB;

    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            // calculate rpart for A
            rPartA = basis.rPart(rVals, m, n);

            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                a = indices(m, n, n_z);
                // calculate zpart for A and funcA
                zPartA = basis.zPart(zVals, n_z);
                arma::mat funcA = rPartA * zPartA.t();
                // add to result the computed value on the diagonal at index a
                result += funcA % funcA * rho(a, a);

                for (int np = 0; np < n; np++) {
                    // calculate rpart for B
                    rPartB = basis.rPart(rVals, m, np);

                    for (int n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++) {
                        b = indices(m, np, n_zp);
                        // calculate zpart for B and funcB
                        zPartB = basis.zPart(zVals, n_zp);
                        arma::mat funcB = rPartB * zPartB.t();
                        // add to result the computed value twice (because rho is diagonal)
                        result += 2 * funcA % funcB * rho(a, b);
                    } // n_zp
                } // np

                // Do one more pass of n_zp for when np == n.
                for (int n_zp = 0; n_zp < n_z; n_zp++) {
                    b = indices(m, n, n_zp);
                    zPartB = basis.zPart(zVals, n_zp);
                    // We can use rPartA here because np == n, thus rPartA == rPartB
                    arma::mat funcB = rPartA * zPartB.t();
                    result += 2 * funcA % funcB * rho(a, b);
                } // n_zp
            } // n_z
        } // n
    } // m
    return result;
}

