#include "../include/Basis.h"
#include "../include/OneDHOSolution.h"
#include "../include/Poly.h"

#include <cmath>
#include <stdexcept>

namespace
{
constexpr double kPi = 3.141592653589793238462643383279502884;
}

Basis::Basis(double br, double bz, int N, double Q)
    : br_(br)
    , bz_(bz)
    , N_(N)
    , Q_(Q)
{
    if (Q_ <= 0.0)
    {
        throw std::invalid_argument("Q must be positive.");
    }
    if (N_ <= 0)
    {
        throw std::invalid_argument("N must be a positive integer.");
    }

    mMax = compute_mMax(N_, Q_);
    nMax = compute_nMax();
    n_zMax = compute_nzMax(N_, Q_);
}

double Basis::compute_Nu(unsigned int i, double N, double Q)
{
    if (Q <= 0.0)
    {
        throw std::invalid_argument("Q must be positive.");
    }
    if (N <= 0.0)
    {
        throw std::invalid_argument("N must be positive.");
    }

    return (N + 2.0) * std::pow(Q, 2.0 / 3.0) + 0.5 - i * Q;
}

int Basis::compute_mMax(double N, double Q)
{
    const double nu0 = compute_Nu(0, N, Q);

    if (nu0 < 1.0)
    {
        throw std::domain_error("No positive index i satisfies nu(i) >= 1.");
    }

    int i = 1;

    // Find the largest integer i such that nu(i) >= 1 with nu being a decreasing function of i
    while (compute_Nu(i, N, Q) >= 1.0)
    {
        ++i;
    }

    return i - 1;
}

arma::ivec Basis::compute_nMax()
{
    // Allocate a vector with dimension nMax
    arma::ivec value(this->mMax);

    for (int i = 0; i < this->mMax; ++i)
    {
        double resultval = 0.5 * (this->mMax - i - 1.0) + 1.0;

        if (resultval < 0.0)
        {
            throw std::domain_error("Computed nMax element is negative; check input parameters.");
        }

        value(i) = resultval;
    }

    return value;
}

arma::imat Basis::compute_nzMax(double N, double Q)
{
    // Allocate matrix with dimensions (mMax) times (nMax)
    int max_of_nMax = 0;
    for (unsigned int i = 0; i < this->nMax.n_elem; ++i)
    {
        if (this->nMax(i) > max_of_nMax)
        {
            max_of_nMax = this->nMax(i);
        }
    }

    arma::imat value(this->mMax, max_of_nMax, arma::fill::zeros);

    for (int m = 0; m < this->mMax; ++m)
    {
        for (int n = 0; n < this->nMax(m); ++n)
        {
            double arg = m + 2.0 * n + 1.0;

            double nuval = compute_Nu(arg, this->N_, this->Q_);

            if (nuval < 0.0)
            {
                throw std::domain_error("Computed nzMax element is negative; check input parameters.");
            }

            value(m, n) = nuval;
        }
    }

    return value;
}

arma::vec Basis::rPart(const arma::vec & r_perps, int m, int n) const
{
    if (n < 0)
    {
        throw std::invalid_argument("n must be non-negative");
    }

    double factorial1 = std::tgamma(static_cast<double>(n) + 1.0);
    double factorial2 = std::tgamma(static_cast<double>(n) + static_cast<double>(std::abs(m)) + 1.0);

    double coeff1 = (1.0 / (br_ * std::sqrt(kPi))) * std::sqrt(factorial1 / factorial2);

    arma::vec coeff2 = arma::exp(-arma::square(r_perps) / (2.0 * br_ * br_));

    const int abs_m = std::abs(m);
    arma::vec coeff3 = arma::pow(r_perps / br_, abs_m);

    arma::vec r_perpstemps = arma::square(r_perps) / (br_ * br_);

    Poly laguerreHelper;
    laguerreHelper.calcLaguerre(abs_m + 1, n + 1, r_perpstemps);
    arma::vec laguerreVals = laguerreHelper.laguerre(abs_m, n);

    return coeff1 * coeff2 % coeff3 % laguerreVals;
}

arma::vec Basis::zPart(const arma::vec & zValues, int n_z) const
{
	if (n_z < 0)
	{
		throw std::invalid_argument("n_z must be non-negative");
	}

	if (bz_ <= 0.0)
	{
		throw std::invalid_argument("bz_ must be positive");
	}

    double factorial = std::tgamma(static_cast<double>(n_z) + 1.0);

    double coeff = 1.0 / (std::sqrt(bz_) * std::pow(kPi, 0.25) * std::sqrt(std::pow(2.0, n_z) * factorial));

    const arma::vec gaussian = arma::exp(-arma::square(zValues) / (2.0 * bz_ * bz_));
    const arma::vec normalizedZ = zValues / bz_;

    Poly hermiteHelper;
    hermiteHelper.calcHermite(n_z, normalizedZ);
    arma::vec hermiteVals = hermiteHelper.getHermiteRow(n_z);

    return coeff * gaussian % hermiteVals;
}

arma::mat Basis::basisFunc(int m, int n, int n_z, const arma::vec& zVals, const arma::vec& rVals) const
{
    arma::vec radial = rPart(rVals, m, n);
    arma::vec axial = zPart(zVals, n_z);
    return radial * axial.t();
}