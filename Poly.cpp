#include "Poly.h"




void Poly::calcHermit( int n , arma::rowvec zVals ) {

    nombreValeur = zVals.n_elem;
    hermit = arma::mat(n, nombreValeur, arma::fill::zeros);

    hermit.row(0).fill(1);
    if (ordreMax > 1) {
        hermit.row(1) = 2 * zVals;

        for (int ligne = 2 ; ligne < ordreMax + 1 ; ligne++ ){
        hermit.row(ligne) = 2.0 * (zVals % hermit.row(ligne - 1)) - 2.0 * (ligne - 1) * hermit.row(ligne - 2);
        }
    }
}




void Poly::calcLaguerre( int n , int m , arma::rowvec zVals ) {


    // =================Partie complétion slice 0=========================
    nombreValeur = zVals.n_elem;
    laguerre = arma::cube(m , nombreValeur,n, arma::fill::zeros);

    printf("test1");
    arma::mat slice0 = arma::mat(m, nombreValeur, arma::fill::zeros); // voir si y'a pas moyen de remplir cette slice directement de 1
    slice0.fill(1);
    laguerre.slice(0)= slice0;
    printf("test2");

    // ===============================================================


 
    // =================Partie complétion slice 1=========================
   //je fais un vect de 1+m 
    //pui pour chaque je  fais un slice avec ce vecteur ajouter à chaqu, évite double boucle


    arma::vec vecTemporaire(m);
    for (int i =0 ; i<m ; i++)
    {

           vecTemporaire(i) = 1 + i;
    }

    arma::mat slice1 = arma::mat(m, nombreValeur);
    for (int nu = 0 ; nu < nombreValeur ; nu++){
        
        arma::vec tempo(m);
        tempo.fill(nu);
        slice1.row(nu) = vecTemporaire + tempo;
    }
    laguerre.slice(1) = slice1;   
    // =========================== Partie complétion autres Slices ============================

    for (int N = 2 ; N < n ; N++){
        laguerre.slice(N) = Poly::calcSliceN(N , laguerre.slice(N-1) , laguerre.slice(N-1));
    }



    // ===============================================================

    Poly::printMatrix(laguerre.slice(0));
    Poly::printMatrix(laguerre.slice(1));
}




void Poly::printMatrix( arma::mat mat ) {
    int width = mat.n_rows;
    int length = mat.n_cols;

    std::cout << "[";
    for (int i = 0 ; i<length ; i++){
        std::cout << "[";
        for (int j = 0 ; j<width ; j++){ 
            std::cout << j << " ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
}

arma::mat Poly::calcSliceN(int n, arma::mat sliceNMoins1, arma::mat sliceNMoins2)
{
    return arma::mat();
}
