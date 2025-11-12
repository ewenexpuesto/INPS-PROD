#ifndef HERMIT_H
#define HERMIT_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

class Hermit
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::rowvec zValues;   
    arma::mat matrixValues; 

    Hermit(int ordreN, const arma::rowvec& zValues_);

    void CalcMatrix();

    arma::rowvec GetHermit();
};

#endif
