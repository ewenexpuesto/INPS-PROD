#ifndef HERMIT_H
#define HERMIT_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Helper for generating Hermite polynomial matrices.
 */
class Hermit
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::rowvec zValues;   
    arma::mat matrixValues; 

    /**
     * @brief Construct the Hermite table for the provided abscissas.
     * @param ordreN highest polynomial order to compute.
     * @param zValues_ row vector containing the sampling points.
     */
    Hermit(int ordreN, const arma::rowvec& zValues_);

    /**
     * @brief Populate matrixValues with Hermite evaluations up to ordreMax.
     */
    void CalcMatrix();

    /**
     * @brief Return the final Hermite row (order ordreMax).
     */
    arma::rowvec GetHermit();
};

#endif
