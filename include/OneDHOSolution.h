#ifndef ONEDHOSOLUTION_H
#define ONEDHOSOLUTION_H

#include "../armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Closed-form solution utilities for the 1D harmonic oscillator.
 */
class OneDHOSolution
{
public:
    double omega = 1.0;
    double hbarre = 1.0;
    double m = 1.0;

    /**
     * @brief Evaluate the nth Hermite polynomial on the provided grid.
     * @param n target polynomial order.
     * @param vecteurAbscisses sample locations (row vector).
     * @param nbEchantillons number of samples to consider.
     * @return Row vector containing the evaluated samples.
     */
    arma::rowvec OneDHOSolutionCalc(int n , arma::rowvec vecteurAbscisses , int nbEchantillons);

    /**
     * @brief Print a solution vector to stdout (debug helper).
     */
    static void PrintOneDHOSolution(const arma::vec vec);
};


#endif //BARBIIE_ONEDHOSOLUTION_H