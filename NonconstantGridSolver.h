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
    Parameters params;
    arma::cx_mat Id;
    int n;
public:
    /*!
     * \brief calculates \f$ \mathbf{Q} \f$ matrix.
     */
    arma::cx_mat Q(double x, double E);

    /*!
     * \brief calculates \f$ \mathbf{T}\f$ matrix for given x (uses interpolation from Parameters).
     */
    arma::cx_mat T(double x, double E, double dx);

    /*!
     * \brief Calculates \f$ \mathbf{U}\f$ matrix for given x (uses interpolation from Parameters).
     */
    arma::cx_mat U(double x, double E, double dx);

    /*!
     * \brief Modifies closed channels elements.
     */
    void modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E);

    /*!
   * \brief Calculates \f$ \mathbf{E}^+ (x, E)\f$ matrix.
   */
    arma::cx_mat calculateEP(double x, double E);

    /*!
* \brief Calculates \f$ \mathbf{E}^- (x, E)\f$ matrix.
*/
    arma::cx_mat calculateEM(double x, double E);
    /*!
     * Calculates the grid points.
     */
    std::vector<double> generateGrid(double E);

    /*!
     * \brief Iterates Numerov algorithm forward up to \f$N-1\f$
     */
    arma::cx_mat fwdIteration(const arma::cx_mat &B, double E, std::vector<double> grid);

    /*!
     * \brief Calculates \f$ \mathbf{S}\f$ matrix given \f$\mathbf{R}_{N-1}
     */
    arma::cx_mat calculateS(const arma::cx_mat R_N, double E, std::vector<double> grid);

    /*!
     * \brief Saves \f$ \mathbf{S}\f$ matrix (Im and Re part separately)
     */
    void saveS(const arma::cx_mat &S, const std::string S_type, const double E);

    void setParameters(const Parameters &parameters) {
        params = parameters;
    }


    NonconstantGridSolver() = default;

    NonconstantGridSolver(const NonconstantGridSolver &) = default;

    ~NonconstantGridSolver() = default;

    /*!
     * \brief Constructor.
     */
    explicit NonconstantGridSolver(const Parameters &params) : params(params) {
        n = params.getNChannels();
        Id = arma::cx_mat(n, n, arma::fill::eye);
    };

    /*!
     * \brief Performs Numerov calculations for a given set of parameters.
     */
    void solveForEnergies(std::string directory);

};


#endif //NOMICROMOTION_NONCONSTANTGRIDSOLVER_H
