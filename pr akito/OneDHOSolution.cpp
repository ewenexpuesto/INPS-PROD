#include "OneDHOSolution.h"
#include "Hermit.h"
#include <cmath>
#include <cstdio>

// Fonction pour calculer la factorielle
double factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

arma::rowvec OneDHOSolution::OneDHOSolutionCalc(int n , arma::rowvec vecteurAbscises , int nbEchantillons)
{
    
    //Calcul de la constante
    double facteurConstant = 1.0 / sqrt(pow(2.0, n) * factorial(n)) * pow((m * omega) / (M_PI * hbarre), 0.25);    
    
    //Vecteurs des coefficients pour chaque z
    arma::vec vec(nbEchantillons);
    arma::rowvec vecteurCoefExp(nbEchantillons);
    for (int i = 0; i < nbEchantillons; i++) vecteurCoefExp(i) = exp(-(omega * m * pow(vecteurAbscises(i), 2)) / (2.0 * hbarre));
    arma::rowvec coefs = facteurConstant * vecteurCoefExp;
        

    Hermit hermit(n , vecteurAbscises);
    arma::rowvec hermiteVals = hermit.GetHermit();
    arma::rowvec vecResults = coefs % hermiteVals;; 
    return vecResults;

}

void OneDHOSolution::PrintOneDHOSolution(const arma::vec vec)
{
    std::cout << "[";
    for (double v : vec)
    {
        std::cout << v << " ";
    }
    std::cout << "]" << std::endl;
}

