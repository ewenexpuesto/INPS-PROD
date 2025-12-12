#include "Poly.h"

#include <stdexcept>


void Poly::calcHermite(int n, const arma::vec & zVals)
{
    if (n < 0)
    {
        throw std::invalid_argument("Hermite order must be non-negative");
    }

    ordreMax = n;
    nombreValeur = zVals.n_elem;
    arma::rowvec zRow = zVals.t();

    // Allocate the Hermite table and seed the constant row
    hermiteTable = arma::mat(n + 1, nombreValeur, arma::fill::zeros);
    hermiteTable.row(0).ones();

    if (n == 0)
    {
        return;
    }

    // Base recurrence row for n = 1
    hermiteTable.row(1) = 2.0 * zRow;

    // Fill higher orders via the standard Hermite recurrence H_n = 2 x H_{n-1} - 2(n-1) H_{n-2}
    for (int ligne = 2; ligne <= n; ++ligne)
    {
        hermiteTable.row(ligne) = 2.0 * (zRow % hermiteTable.row(ligne - 1)) - 2.0 * (ligne - 1) * hermiteTable.row(ligne - 2);
    }
}


void Poly::calcLaguerre(int mCount, int nCount, const arma::vec& zVals)
{
    nombreValeur = zVals.n_elem;

    if (mCount <= 0 || nCount <= 0)
    {
        laguerreTable.reset();
        return;
    }

    // Allocate the cube: rows=m alphas, cols=n orders, slices=x samples
    laguerreTable = arma::cube(mCount, nCount, nombreValeur, arma::fill::zeros);

    for (int alpha = 0; alpha < mCount; ++alpha)
    {
        const double alpha_d = static_cast<double>(alpha);

        for (int idx = 0; idx < nombreValeur; ++idx)
        {
            const double x = zVals(idx);

            // L_0^alpha(x) = 1 for all x
            laguerreTable(alpha, 0, idx) = 1.0;

            if (nCount == 1)
            {
                continue;
            }

            // L_1^alpha(x) = 1 + alpha - x
            laguerreTable(alpha, 1, idx) = 1.0 + alpha_d - x;

            // Higher orders use the generalized Laguerre recurrence
            for (int order = 2; order < nCount; ++order)
            {
                const double order_d = static_cast<double>(order);
                const double coeffA = 2.0 * order_d - 1.0 + alpha_d - x;
                const double coeffB = order_d - 1.0 + alpha_d;
                laguerreTable(alpha, order, idx) = (coeffA * laguerreTable(alpha, order - 1, idx) - coeffB * laguerreTable(alpha, order - 2, idx)) / order_d;
            }
        }
    }
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
    arma::mat result(sliceNMoins1.n_rows, sliceNMoins1.n_cols, arma::fill::zeros);
    
    result = sliceNMoins1;

    return result;
}

arma::vec Poly::getHermiteRow(int order) const
{
    if (hermiteTable.is_empty())
    {
        throw std::logic_error("Hermite table is empty; call calcHermite first");
    }

    return hermiteTable.row(order).t();
}

arma::vec Poly::hermite(int order) const
{
    return getHermiteRow(order);
}

arma::vec Poly::laguerre(int mOrder, int nOrder) const
{
    // Extract the tube (all x values) for the requested (m,n) pair
    return laguerreTable.tube(mOrder, nOrder);
}
