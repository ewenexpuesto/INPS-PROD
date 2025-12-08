
#include <armadillo>
#include "basis.h"

// Solver class to calculate the nuclear density with all different level of optimisation 

class Solver {
  public:
    //Default constructor
    Solver(Basis);
    void reset();
    //naive implementation
    arma::mat version0(arma::Col<double>, arma::Col<double>);

    //same as version0 but that removes a loop
    arma::mat version1(arma::Col<double>, arma::Col<double>);

    //same as version1 but that take into account the fact that rho is symmetric
    arma::mat version2(arma::Col<double>, arma::Col<double>);

    //same as version2 but calculating some values higher in the nested loops
    arma::mat version3(arma::Col<double>, arma::Col<double>);

    //same as version3 but removing all of the conditional tests except the ones in function calls
    arma::mat version4(arma::Col<double>, arma::Col<double>);

    //same as version4 but not calculating moving calls to rPart higher up in the nested loop
    arma::mat version5(arma::Col<double>, arma::Col<double>);

    /// The basis used
    Basis basis;

    /// Matrix to store rho values
    arma::mat rho;

    /// Matrix to store indices linked to unique ( \a m, \a n, \a nz ) quantum numbers
    arma::icube indices;

};

#endif


