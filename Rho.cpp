#include "Rho.h"

arma::mat Rho::compute(const Basis& basis,
					   const arma::mat& rhoMatrix,
					   const arma::vec& rVals,
					   const arma::vec& zVals)
{
	arma::mat result = arma::zeros(rVals.n_elem, zVals.n_elem);

	arma::uword idxA = 0;
	for (int m = 0; m < basis.mMax; ++m)
	{
		for (int n = 0; n < basis.nMax(m); ++n)
		{
			for (int n_z = 0; n_z < basis.n_zMax(m, n); ++n_z)
			{
				arma::mat funcA = basis.basisFunc(m, n, n_z, zVals, rVals);

				arma::uword idxB = 0;
				for (int mp = 0; mp < basis.mMax; ++mp)
				{
					for (int np = 0; np < basis.nMax(mp); ++np)
					{
						for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); ++n_zp)
						{
							arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
							double coeff = rhoMatrix(idxA, idxB);
							result += funcA % funcB * coeff;
							++idxB;
						}
					}
				}

				++idxA;
			}
		}
	}

	return result;
}
