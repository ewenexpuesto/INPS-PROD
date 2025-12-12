#ifndef SOLVER_H
#define SOLVER_H

#include "Basis.h"

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

// Solver class to calculate the nuclear density with all different level of optimisation 

class Solver {
  public:
    //Default constructor
    Solver(const Basis&);
    void reset();
    //naive implementation
    arma::mat version0(const arma::vec&, const arma::vec&);

    //same as version0 but that removes a loop
    arma::mat version1(const arma::vec&, const arma::vec&);

    //same as version1 but that take into account the fact that rho is symmetric
    arma::mat version2(const arma::vec&, const arma::vec&);

    //same as version2 but calculating some values higher in the nested loops
    arma::mat version3(const arma::vec&, const arma::vec&);

    //same as version3 but removing all of the conditional tests except the ones in function calls
    arma::mat version4(const arma::vec&, const arma::vec&);

    //same as version4 but not calculating moving calls to rPart higher up in the nested loop
    arma::mat version5(const arma::vec&, const arma::vec&);

    /// The Basis used
    Basis basis;

    /// Matrix to store rho values
    arma::mat rho;

    /// Matrix to store indices linked to unique ( \a m, \a n, \a nz ) quantum numbers
    arma::icube indices;

};

#endif


