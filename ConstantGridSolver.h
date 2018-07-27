//
// Created by martas on 07.12.17.
//

#ifndef NOMICROMOTION_CONSTANTGRIDSOLVER_H
#define NOMICROMOTION_CONSTANTGRIDSOLVER_H

/*!
 * \file
 * \brief Definition of Solver class
 *
 * This file contains a definition of Solver class, performing the calculations for
 * a given set of parameters (Parameters object).
 */

#include "armadillo"
#include "Parameters.h"
#include <sstream>
#include <ostream>
#include <fstream>
#include <cstring>
#include <math.h>
#include <gtest/gtest_prod.h>



class ConstantGridSolver{


    /*!
     * \brief Set of parameters for which the calculations are performed.
     */
    Parameters params;
public:
    /*!
     * \brief Calculates T matrix for given index
     */
    arma::cx_mat calculateT(int j, double E) const;
    /*!
     * \brief Calculates U matrix
     */
    arma::cx_mat calculateU(int j, double E);

    /*!
     * \brief Calculates EP matrices
     */
    arma::cx_mat calculateEP(int ind_x, double E);
    /*!
     * \brief Calculates EM matrices
     */
    arma::cx_mat calculateEM(int ind_x, double E);
    /*!
     * \brief Modifies closed channels elements
     */
    void modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E);
    /*!
     * \brief Iterates Numerov algorithm forward up to N
     */
    arma::cx_mat fwdIteration(const arma::cx_mat &B, double E);
    /*!
     * \brief Calculates S matrix given RN
     */
    arma::cx_mat calculateS(const arma::cx_mat R_N, double E);
    /*!
     * \brief Saves S matrix (Im and Re part separately)
     */
    void saveS(const arma::cx_mat &S, const std::string S_type, const double E, const std::string directory);

    void setParameters(const Parameters &parameters){
        params = parameters;
    }



    ConstantGridSolver() = default;
   // ConstantGridSolver(const ConstantGridSolver&) = default;
    ~ConstantGridSolver() = default;

    /*!
     * \brief Constructor.
     */
     explicit ConstantGridSolver(const Parameters& params);
    /*!
     * \brief Performs Numerov calculations for a given set of parameters.
     */
    void solveForEnergies(std::string directory);

};


#endif //NOMICROMOTION_SOLVER_H
