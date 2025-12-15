#include "../include/Hermit.h"

Hermit::Hermit(int ordreN, const arma::rowvec& zValues_) {
    ordreMax = ordreN;
    zValues = zValues_;
    nombreValeur = zValues_.n_elem;
    matrixValues = arma::mat(ordreMax +1, nombreValeur, arma::fill::zeros);
}

void Hermit::CalcMatrix() {

    matrixValues.row(0).fill(1);

    if (ordreMax > 1) {

        matrixValues.row(1) = 2 * zValues;

        for (int ligne = 2 ; ligne < ordreMax + 1 ; ligne++ ){

        matrixValues.row(ligne) = 2.0 * (zValues % matrixValues.row(ligne - 1)) - 2.0 * (ligne - 1) * matrixValues.row(ligne - 2);
        }
    }

}
arma::rowvec Hermit::GetHermit() {
    CalcMatrix();
    return matrixValues.row(ordreMax);
}
