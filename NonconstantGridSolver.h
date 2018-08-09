//
// Created by martas on 22.06.18.
//

#ifndef NOMICROMOTION_NONCONSTANTGRIDSOLVER_H
#define NOMICROMOTION_NONCONSTANTGRIDSOLVER_H

#include "armadillo"
#include "Parameters.h"
#include <sstream>
#include <ostream>
#include <fstream>
#include <cstring>
#include <math.h>

class NonconstantGridSolver {
    /*!
     * \brief Set of parameters for which the calculations are performed.
     */
    const Parameters &params;
    arma::cx_mat Id;
    int n;
    std::vector<double> gridPoints;
    /*!
     * \brief calculates T matrix for given x (uses interpolation from Parameters)
     */
    arma::cx_mat calculateT(double x, double E, double dx);
    /*!
     * \brief Calculates U matrix for given x (uses interpolation from Parameters)
     */
    arma::cx_mat calculateU(double x, double E, double dx);
    /*!
     * \brief Calculates F(alpha) based on interpolation formula, needed for halving the grid. For forward iteration.
     */
    arma::cx_mat
            calculateFAlpha(double x, double dx, double E, const arma::cx_mat &F_i, const arma::cx_mat &F_i1,
                            const arma::cx_mat &U_i);
    /*!
     * \brief Modifies closed channels elements
     */
    void modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E);
    /*!
     * \brief Calculates EM matrices
     */
    auto calculateEMP(double x0, double x1, double E);

    /*!
     * \brief Iterates Numerov algorithm forward up to N
     */
    arma::cx_mat fwdIteration(const arma::cx_mat &B, double E);
    /*!
     * \brief Calculates S matrix given RN
     */
    arma::cx_mat calculateS(const arma::cx_mat R_N, double E, double dx);
    /*!
     * \brief Saves S matrix (Im and Re part separately)
     */
    void saveS(const arma::cx_mat &S, const std::string S_type, const double E);
    /*!
     * \brief Checks what is the required grid in point x and how it compares to the current value.
     * Returns halved/double/the same value.
     */
    double changeDxIfNeeded(double xBef, double dxBef, double E);


public:
    NonconstantGridSolver() = delete;
    NonconstantGridSolver(const NonconstantGridSolver&) = delete;
    ~NonconstantGridSolver() = default;
    /*!
     * \brief Constructor.
     */
    explicit NonconstantGridSolver(const Parameters& params):params(params){
        n = params.getNChannels();
        std::cout << n << '\n';
        Id = arma::cx_mat(n, n, arma::fill::eye);
        Id.print();
    };
    /*!
     * \brief Performs Numerov calculations for a given set of parameters.
     */
    void solveForEnergies();

};


#endif //NOMICROMOTION_NONCONSTANTGRIDSOLVER_H
