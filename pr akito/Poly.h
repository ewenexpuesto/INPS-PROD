#ifndef POLY_H
#define POLY_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

class Hermit
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::mat hermit; 
    arma::mat laguerre;
    
    void calcHermit(int n , arma::rowvec zVals);
    void calcLaguerre(int n , int m , arma::rowvec zVals);


};

#endif
