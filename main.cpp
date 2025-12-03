#include "Poly.h"
#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"



int main(){
    Poly poly;
    printf("test");
    arma::rowvec v = arma::randu<arma::rowvec>(10); // 10 réels uniformes entre 0 et 1    Poly poly;
    poly.calcLaguerre(5 , 6 , v);
    return (0);
}