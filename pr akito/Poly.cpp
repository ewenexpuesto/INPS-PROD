#include "Hermit.h"



void Poly::calcHermit( int n , arma::rowvec zVals ) {

    nombreValeur = zVals.n_elem;
    hermit = arma::mat(n, nombreValeur, arma::fill::zeros);


    hermit.row(0).fill(1);
    if (ordreMax > 1) {
        hermit.row(1) = 2 * zValues;

        for (int ligne = 2 ; ligne < ordreMax + 1 ; ligne++ ){
        hermit.row(ligne) = 2.0 * (zValues % hermit.row(ligne - 1)) - 2.0 * (ligne - 1) * hermit.row(ligne - 2);
        }
    }


}
arma::rowvec Hermit::GetHermit() {

}
