//
// Created by martas on 07.12.17.
//

#ifndef NOMICROMOTION_CONSTANTGRIDSOLVER_H
#define NOMICROMOTION_CONSTANTGRIDSOLVER_H

/*!
 * \file
 * \brief Definition of ConstantGridSolver class
 *
 * This file contains a definition of ConstantGridSolver class, performing the calculations for
 * a given set of parameters (Parameters object).
 */

#include "armadillo"
#include "Parameters.h"
#include <sstream>
#include <ostream>
#include <fstream>
#include <cstring>
#include <math.h>


class ConstantGridSolver {
    /*!
     * \brief Set of parameters for which the calculations are performed.
     */
    Parameters params;
public:
    /*!
     * \brief Calculates \f$ \mathbf{T}(x_j, E) \f$ matrix.
     */
    arma::cx_mat calculateT(int j, double E) const;

    /*!
     * \brief Calculates \f$ \mathbf{U}(x_j, E) \f$ matrix.
     */
    arma::cx_mat calculateU(int j, double E);

    /*!
     * \brief Calculates \f$ \mathbf{E}^+ (x_j, E)\f$ matrix.
     */
    arma::cx_mat calculateEP(int j, double E);

    /*!
     * \brief Calculates \f$ \mathbf{E}^-(x_j, E)\f$ matrix.
     */
    arma::cx_mat calculateEM(int j, double E);

    /*!
     * \brief Modifies closed channels elements.
     */
    void modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E);

    /*!
     * \brief Iterates Numerov algorithm forward up to \f$ N-1 \f$ and returns \f$ \mathbf{R}_{N-1} \f$ matrix for a given energy.
     */
    arma::cx_mat fwdIteration(const arma::cx_mat &B, double E);

    /*!
     * \brief Calculates \f$ \mathbf{S} \f$ matrix for given \f$ \mathbf{R}_{N-1} \f$
     */
    arma::cx_mat calculateS(const arma::cx_mat R_N, double E);

    /*!
     * \brief Saves \f$ \mathbf{S} \f$ matrix (Im and Re part separately).
     */
    void saveS(const arma::cx_mat &S, const int E, const std::string directory);

    void setParameters(const Parameters &parameters) {
        params = parameters;
    }


    ConstantGridSolver() = default;

    // ConstantGridSolver(const ConstantGridSolver&) = default;
    ~ConstantGridSolver() = default;

    /*!
     * \brief Constructor.
     */
    explicit ConstantGridSolver(const Parameters &params);

    /*!
     * \brief Performs Numerov calculations for a given set of parameters for all energies.
     */
    void solveForEnergies(std::string directory);

};


#endif //NOMICROMOTION_SOLVER_H
