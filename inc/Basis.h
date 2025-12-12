#ifndef BASIS_H
#define BASIS_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Encapsulates the truncated cylindrical harmonic oscillator basis.
 *
 * Provides helper routines for evaluating radial (rPart) and axial (zPart)
 * components along with precomputed truncation tables derived from the input
 * quantum numbers.
 */
class Basis
{
public:
    /**
     * @brief Construct the basis descriptors and truncation tables.
     * @param br radial oscillator length scale.
     * @param bz axial oscillator length scale.
     * @param N  principal truncation parameter.
     * @param Q  anisotropy ratio, must be strictly positive.
     */
    Basis(double br, double bz, int N, double Q);

    /**
     * @brief Evaluate the radial portion of the basis at the given radii.
     * @param r_perps vector of radial coordinates.
     * @param m azimuthal quantum number.
     * @param n radial quantum number.
     * @return Normalized radial wavefunction samples.
     */
    arma::vec rPart(const arma::vec& r_perps, int m, int n) const;

    /**
     * @brief Evaluate the axial portion of the basis at the given z values.
     * @param zValues vector of axial coordinates.
     * @param n_z axial quantum number (non-negative).
     * @return Normalized axial wavefunction sapresmples.
     */
    arma::vec zPart(const arma::vec& zValues, int n_z) const;

    /**
     * @brief Build the full basis function ψ_{m,n,nz}(r,z).
     * @param m Azimuthal quantum number.
     * @param n Radial quantum number.
     * @param n_z Axial quantum number.
     * @param zVals z-axis sampling points.
     * @param rVals radial sampling points.
     * @return Outer product R(r) * Z(z)^T.
     */
    arma::mat basisFunc(int m, int n, int n_z, const arma::vec& zVals, const arma::vec& rVals) const;

    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;

private:
    double br_;
    double bz_;
    int N_;
    double Q_;

    /**
     * @brief Helper computing the nu(i) truncation boundary.
     */
    static double compute_Nu(unsigned int i, double N, double Q);

    /**
     * @brief Determine the largest admissible azimuthal index.
     */
    int compute_mMax(double N, double Q);

    /**
     * @brief Build the per-m radial truncation vector.
     */
    arma::ivec compute_nMax();

    /**
     * @brief Build the per-(m,n) axial truncation matrix.
     */
    arma::imat compute_nzMax(double N, double Q);
};

#endif // BASIS_H
