#ifndef BARBIIE_ONEDHOSOLUTION_H
#define BARBIIE_ONEDHOSOLUTION_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"


class OneDHOSolution
{
public:
    double omega = 1.0;
    double hbarre = 1.0;
    double m = 1.0;

    
    //Calcule le polynome d'Hermite d'ordre n, en z 
    arma::rowvec OneDHOSolutionCalc(int n , arma::rowvec vecteurAbscisses , int nbEchantillons);

    static void PrintOneDHOSolution(const arma::vec vec);
};


#endif //BARBIIE_ONEDHOSOLUTION_H