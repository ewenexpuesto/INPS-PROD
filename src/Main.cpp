// #include "Poly.h"
#include "Basis.h"
#include "Poly.h"
#include "Solver.h"
#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#define TS_ASSERT_EQUALS(lhs, rhs) assert((lhs) == (rhs))
#define TS_ASSERT(expr) assert((expr))
#define TS_ASSERT_DELTA(lhs, rhs, delta) assert(std::abs((lhs) - (rhs)) <= (delta))

std::string cubeToDf3(const arma::cube& m)
{
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    int nx = static_cast<int>(m.n_rows);
    int ny = static_cast<int>(m.n_cols);
    int nz = static_cast<int>(m.n_slices);
    ss.put(static_cast<char>(nx >> 8));
    ss.put(static_cast<char>(nx & 0xff));
    ss.put(static_cast<char>(ny >> 8));
    ss.put(static_cast<char>(ny & 0xff));
    ss.put(static_cast<char>(nz >> 8));
    ss.put(static_cast<char>(nz & 0xff));
    double theMin = 0.0;
    double theMax = m.max();
    double denom = theMax - theMin;
    if (denom == 0.0)
    {
        denom = 1.0;
    }
    for (arma::uword k = 0; k < m.n_slices; ++k)
    {
        for (arma::uword j = 0; j < m.n_cols; ++j)
        {
            for (arma::uword i = 0; i < m.n_rows; ++i)
            {
                unsigned int v = static_cast<unsigned int>(
                    255.0 * (std::fabs(m(i, j, k)) - theMin) / denom);
                ss.put(static_cast<char>(v));
            }
        }
    }
    return ss.str();
}

void saveCubeAsDf3(const arma::cube& cube, const std::string& filename)
{
    std::string df3Data = cubeToDf3(cube);
    std::ofstream out(filename, std::ios::binary);
    if (!out)
    {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    out.write(df3Data.data(), static_cast<std::streamsize>(df3Data.size()));
    std::cout << filename << " written (" << df3Data.size() << " bytes)" << std::endl;
}

int main()
{
    // Poly poly;
    // printf("test");
    // arma::rowvec v = arma::randu<arma::rowvec>(10); // 10 réels uniformes entre 0 et 1    Poly poly;
    // poly.calcLaguerre(5 , 6 , v);
    // return (0);

    // Mandatory test #00 - Hermite and Laguerre polynomials
    Poly poly;
    arma::vec zVals, calcVals, targetVals;
    zVals = {-3.1, -2.3, -1.0, -0.3, 0.1, 4.3, 9.2, 13.7};
    poly.calcHermite(6, zVals); // compute Hermite polynomials for n in {0 ... 5}
    calcVals   = poly.hermite(4); // n = 4
    targetVals = {  1.02835360e+03,  2.05825600e+02, -2.00000000e+01,  7.80960000e+00,
                    1.15216000e+01,  4.59456160e+03,  1.10572154e+05,  5.54643458e+05};
    TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);
    calcVals   = poly.hermite(5); // n = 5
    targetVals = { -4.76676832e+03, -3.88909760e+02,  8.00000000e+00, -3.17577600e+01,
                    1.18403200e+01,  3.48375818e+04,  1.98557479e+06,  1.50339793e+07};
    TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);
    zVals = {0.1, 0.3, 1.2, 1.8, 2.0, 2.5, 7.1, 11.1};
    poly.calcLaguerre(6, 4, zVals); // compute generalized Laguerre polynomials for m in {0 ... 5} and n in {0 ... 3}
    calcVals   = poly.laguerre(4, 2); // m = 4, n = 2
    targetVals = {  14.405,  13.245,  8.52 ,  5.82 ,  5.,  3.125,  -2.395,  10.005};
    TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);
    calcVals   = poly.laguerre(5, 3); // m = 5, n = 3
    targetVals = { 53.23983333,  47.95550000,  27.87200000,  17.5880,
                14.66666667,   8.39583333,  -0.81183333,  10.1015};
    TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);

    printf("Hermite and Laguerre polynomials tests passed!\n");

    // Mandatory test #01 - Basis truncation
    //     br = 1.935801664793151, bz = 2.829683956491218, N = 14, Q = 1.3
    Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);
    TS_ASSERT_EQUALS(basis.mMax, 14);
    arma::ivec nMax = {7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1};
    TS_ASSERT((basis.nMax - nMax).is_zero());
    arma::imat n_zMax = {{18, 15, 13, 10, 7, 5, 2}, 
                        {16, 14, 11,  9, 6, 3, 1}, 
                        {15, 13, 10,  7, 5, 2, 0}, 
                        {14, 11,  9,  6, 3, 1, 0}, 
                        {13, 10,  7,  5, 2, 0, 0}, 
                        {11,  9,  6,  3, 1, 0, 0}, 
                        {10,  7,  5,  2, 0, 0, 0}, 
                        { 9,  6,  3,  1, 0, 0, 0}, 
                        { 7,  5,  2,  0, 0, 0, 0}, 
                        { 6,  3,  1,  0, 0, 0, 0}, 
                        { 5,  2,  0,  0, 0, 0, 0}, 
                        { 3,  1,  0,  0, 0, 0, 0}, 
                        { 2,  0,  0,  0, 0, 0, 0}, 
                        { 1,  0,  0,  0, 0, 0, 0}};
    // check if matrices are equal
    TS_ASSERT((basis.n_zMax - n_zMax).is_zero());

    printf("Basis truncation tests passed!\n");

    // Mandatory test #02 - Basis r-functions
    //     br = 1.935801664793151, bz = 2.829683956491218, N = 14, Q = 1.3
    arma::vec r = {3.1, 2.3, 1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
    arma::vec res00 = { 8.08521235111303e-02,
                        1.43887615825118e-01,
                        2.55045100912706e-01,
                        2.91450097294984e-01,
                        2.91061479407116e-01,
                        2.47240792330589e-02,
                        3.63004153921473e-06,
                        3.87659726026123e-12};
    TS_ASSERT_DELTA(arma::norm(basis.rPart(r, 0, 0) - res00), 0.0, 1e-15);
    arma::vec res82 = { 5.87858442372438e-02,
                        1.35240488413384e-02,
                        4.06810074575519e-05,
                        0.00000000000000e+00,
                        4.92817669085478e-13,
                        8.52011998934850e-02,
                        5.20525909328609e-02,
                        1.44615166152252e-05};
    TS_ASSERT_DELTA(arma::norm(basis.rPart(r, 8, 2) - res82), 0.0, 1e-15);

    printf("Basis r-function tests passed!\n");

    // Mandatory test #03 - Basis z-functions
    //     br = 1.935801664793151, bz = 2.829683956491218, N = 14, Q = 1.3
    arma::vec z = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
    arma::vec res00bis = { 7.64546544834383e-04,
                        5.44886272162148e-03,
                        4.19492564268520e-01,
                        4.46522724110539e-01,
                        4.46243982300708e-01,
                        1.40736821086932e-01,
                        2.26186220733178e-03,
                        3.62929640195959e-06};
    TS_ASSERT_DELTA(arma::norm(basis.zPart(z, 0) - res00bis), 0.0, 1e-15);
    arma::vec res15 = {-9.48674551049192e-02,
                    -1.40338701953237e-03,
                        1.85620628040096e-01,
                    -0.00000000000000e+00,
                    -3.93028470685214e-02,
                    -1.79526868763440e-01,
                        2.15604096600475e-01,
                        2.44977220882127e-01};
    TS_ASSERT_DELTA(arma::norm(basis.zPart(z, 15) - res15), 0.0, 1e-15);

    printf("Basis z-function tests passed!\n");

    // Test #04 - Basis size
    int basisVectorCount = 0;
    for (int m = 0; m < basis.mMax; ++m)
    {
        for (int n = 0; n < basis.nMax(m); ++n)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); ++n_z)
            {
                std::cout << "Basis vector " << basisVectorCount
                          << ": m=" << m
                          << " n=" << n
                          << " n_z=" << n_z << std::endl;
                ++basisVectorCount;
            }
        }
    }
    TS_ASSERT_EQUALS(basisVectorCount, 374);
    printf("Basis contains %d vectors as expected.\n", basisVectorCount);
    
    printf("All tests passed!\n");

    // Test #05 - Solver versions (simple benchmarks)
    Solver solver(basis);
    arma::vec visualizationX = arma::linspace<arma::vec>(-10.0, 10.0, 32);
    arma::vec visualizationY = arma::linspace<arma::vec>(-10.0, 10.0, 32);
    arma::vec visualizationZ = arma::linspace<arma::vec>(-20.0, 20.0, 64);
    arma::vec visualizationR(visualizationX.n_elem);
    for (arma::uword i = 0; i < visualizationR.n_elem; ++i)
    {
        double x = visualizationX(i);
        double y = visualizationY(i);
        visualizationR(i) = std::sqrt(x * x + y * y);
    }

    arma::wall_clock wclock;

    wclock.tic();
    arma::mat density0 = solver.version0(visualizationR, visualizationZ);
    double elapsed0 = wclock.toc();
    std::cout << "Solver::version0 -> " << elapsed0 << " s; density size: "
              << density0.n_rows << "x" << density0.n_cols << std::endl;

    wclock.tic();
    arma::mat density1 = solver.version1(visualizationR, visualizationZ);
    double elapsed1 = wclock.toc();
    std::cout << "Solver::version1 -> " << elapsed1 << " s; density size: "
              << density1.n_rows << "x" << density1.n_cols << std::endl;

    wclock.tic();
    arma::mat density2 = solver.version2(visualizationR, visualizationZ);
    double elapsed2 = wclock.toc();
    std::cout << "Solver::version2 -> " << elapsed2 << " s; density size: "
              << density2.n_rows << "x" << density2.n_cols << std::endl;

    wclock.tic();
    arma::mat density3 = solver.version3(visualizationR, visualizationZ);
    double elapsed3 = wclock.toc();
    std::cout << "Solver::version3 -> " << elapsed3 << " s; density size: "
              << density3.n_rows << "x" << density3.n_cols << std::endl;

    wclock.tic();
    arma::mat density4 = solver.version4(visualizationR, visualizationZ);
    double elapsed4 = wclock.toc();
    std::cout << "Solver::version4 -> " << elapsed4 << " s; density size: "
              << density4.n_rows << "x" << density4.n_cols << std::endl;

    wclock.tic();
    arma::mat density5 = solver.version5(visualizationR, visualizationZ);
    double elapsed5 = wclock.toc();
    std::cout << "Solver::version5 -> " << elapsed5 << " s; density size: "
              << density5.n_rows << "x" << density5.n_cols << std::endl;

    // Test #06 - Full 3D density calculation (rho(X,Y,Z) from rho(R,Z))
    std::cout << "Test #06 - Building density cubes from rho(R, Z) grids" << std::endl;

    arma::umat radialRowIndex(visualizationX.n_elem, visualizationY.n_elem, arma::fill::zeros);
    for (arma::uword ix = 0; ix < visualizationX.n_elem; ++ix)
    {
        for (arma::uword iy = 0; iy < visualizationY.n_elem; ++iy)
        {
            double radius = std::sqrt(visualizationX(ix) * visualizationX(ix) + visualizationY(iy) * visualizationY(iy));
            arma::uword closestIdx = 0;
            double bestDiff = std::numeric_limits<double>::max();
            for (arma::uword ir = 0; ir < visualizationR.n_elem; ++ir)
            {
                double diff = std::abs(visualizationR(ir) - radius);
                if (diff < bestDiff)
                {
                    bestDiff = diff;
                    closestIdx = ir;
                }
            }
            radialRowIndex(ix, iy) = closestIdx;
        }
    }

    auto buildCubeFromMatrix = [&](const arma::mat& rhoRZ)
    {
        arma::cube cube(visualizationX.n_elem, visualizationY.n_elem, visualizationZ.n_elem, arma::fill::zeros);
        for (arma::uword ix = 0; ix < visualizationX.n_elem; ++ix)
        {
            for (arma::uword iy = 0; iy < visualizationY.n_elem; ++iy)
            {
                arma::uword rIdx = radialRowIndex(ix, iy);
                for (arma::uword iz = 0; iz < visualizationZ.n_elem; ++iz)
                {
                    cube(ix, iy, iz) = rhoRZ(rIdx, iz);
                }
            }
        }
        return cube;
    };

    wclock.tic();
    arma::cube cube0 = buildCubeFromMatrix(density0);
    double cubeElapsed0 = wclock.toc();
    std::cout << "rho(X,Y,Z) version0 -> " << cubeElapsed0 << " s; cube size: "
              << cube0.n_rows << "x" << cube0.n_cols << "x" << cube0.n_slices << std::endl;
    saveCubeAsDf3(cube0, "density_version0.df3");

    wclock.tic();
    arma::cube cube1 = buildCubeFromMatrix(density1);
    double cubeElapsed1 = wclock.toc();
    std::cout << "rho(X,Y,Z) version1 -> " << cubeElapsed1 << " s; cube size: "
              << cube1.n_rows << "x" << cube1.n_cols << "x" << cube1.n_slices << std::endl;
    saveCubeAsDf3(cube1, "density_version1.df3");

    wclock.tic();
    arma::cube cube2 = buildCubeFromMatrix(density2);
    double cubeElapsed2 = wclock.toc();
    std::cout << "rho(X,Y,Z) version2 -> " << cubeElapsed2 << " s; cube size: "
              << cube2.n_rows << "x" << cube2.n_cols << "x" << cube2.n_slices << std::endl;
    saveCubeAsDf3(cube2, "density_version2.df3");

    wclock.tic();
    arma::cube cube3 = buildCubeFromMatrix(density3);
    double cubeElapsed3 = wclock.toc();
    std::cout << "rho(X,Y,Z) version3 -> " << cubeElapsed3 << " s; cube size: "
              << cube3.n_rows << "x" << cube3.n_cols << "x" << cube3.n_slices << std::endl;
    saveCubeAsDf3(cube3, "density_version3.df3");

    wclock.tic();
    arma::cube cube4 = buildCubeFromMatrix(density4);
    double cubeElapsed4 = wclock.toc();
    std::cout << "rho(X,Y,Z) version4 -> " << cubeElapsed4 << " s; cube size: "
              << cube4.n_rows << "x" << cube4.n_cols << "x" << cube4.n_slices << std::endl;
    saveCubeAsDf3(cube4, "density_version4.df3");

    wclock.tic();
    arma::cube cube5 = buildCubeFromMatrix(density5);
    double cubeElapsed5 = wclock.toc();
    std::cout << "rho(X,Y,Z) version5 -> " << cubeElapsed5 << " s; cube size: "
              << cube5.n_rows << "x" << cube5.n_cols << "x" << cube5.n_slices << std::endl;
    saveCubeAsDf3(cube5, "density_version5.df3");
}