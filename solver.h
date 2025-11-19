
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
}
