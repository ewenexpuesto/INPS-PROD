#ifndef POLY_H
#define POLY_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

class Poly
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::mat hermit; 
    arma::cube laguerre;
    
    void calcHermit(int n , arma::rowvec zVals);
    void calcLaguerre(int n , int m , arma::rowvec zVals);
    void static printMatrix(arma::mat mat);
    arma::mat calcSliceN(int n , arma::mat sliceNMoins1 , arma::mat sliceNMoins2); 


};

#endif
