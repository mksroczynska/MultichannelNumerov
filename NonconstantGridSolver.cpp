//
// Created by martas on 22.06.18.
//

#include "NonconstantGridSolver.h"

/*!
 * This method calculates \f$ \mathbf{Q}\f$ matrix for a given point \f$x\f$ and given energy
 * according to the formula:
 * \f{equation}{
 * \mathbf{Q}(x) = -\frac{1}{unit}\left(\mathbf{V}(x) - E\mathbf{I}\right).
 *  \f}
 * @param x
 * @param E - energy
 * @return \f$ \mathbf{Q}(x, E)\f$
 */
arma::cx_mat NonconstantGridSolver::Q(double x, double E) {
    arma::cx_mat V;
    try { V = params.getV(x); }
    catch (std::invalid_argument &ex) {
        std::cout << ex.what();
    }
    arma::cx_mat result;
    try {
        result = -1. / params.getUnit() * (V - E * params.Id());
    }
    catch (...) {
        throw std::runtime_error("Error in Q(x = " + std::to_string(x) + ", E = " + std::to_string(E) + ")");
    }
    return result;
}

/*!
 * This method calculates \f$ \mathbf{T}\f$ matrix for a given point \f$x\f$ and given energy
 * according to the formula:
 * \f{equation}{
 * \mathbf{T}(x) = -\frac{dx}{12}\mathbf{Q}(x).
 *  \f}
 * @param x
 * @param E - energy
 * @param dx - the distance from the next point on the grid
 * @return \f$ \mathbf{T}(x, E)\f$
 */
arma::cx_mat NonconstantGridSolver::T(double x, double E, double dx) {
    return -(dx * dx * Q(x, E)) / 12.;
}

/*!
 * This method calculates \f$ \mathbf{U} \f$ matrix for given \f$x\f$ using the set of parameters provided to the NonconstantGridSolver
 * object according to the following formula:
 * \f{equation}{
 * \mathbf{U}(x) = 12 (\mathbf{I} - \mathbf{T}(x))^{-1} - 10 \mathbf{I}. \f}
 * \param[in] x
 * \param[in] E - energy value
 *  * @param dx - the distance from the next point on the grid
 * \return \f$ \mathbf{U}(x, E) \f$
 * \throws std::invalid_argument if \f$ x \f$ is out of range
*/
arma::cx_mat NonconstantGridSolver::U(double x, double E, double dx) {
    return 12. * arma::inv(Id - T(x, E, dx)) - 10. * Id;
}

void
NonconstantGridSolver::modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E) {
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        if (!params.isOpen(i, E)) {
            j1(i, i) = j1(i, i) / j0(i, i);
            j0(i, i) = 1;
            n1(i, i) = n1(i, i) / n0(i, i);
            n0(i, i) = 1;
        }
    }
}
/*!
 * This method calculates \f$ \mathbf{E}^- \f$ matrix for for a given point \f$x\f$  and given energy.
 * The matrix is diagonal and its elements are calculated the following way:
 * \f{itemize}{
 * \item $ \mathbf{E}^-_{n, n}(x, E) = \exp(-i k x)$ if channel $n$ is open
 * \item $ \mathbf{E}^-_{n, n}(x, E) = \cosh(k x)$ if channel $n$ is closed
 * \item $ \mathbf{E}^-_{n, m}(x, E) = 0$ for $n\neq m$
 * \f}
 * @param x
 * @param E - energy
 * @return \f$ \mathbf{E}^-(x, E)\f$
 * @throws std::invalid_argument if j is wrong
 */
arma::cx_mat NonconstantGridSolver::calculateEM(double x, double E) {
    arma::cx_mat em = params.Id();
    double k_x0;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        try {
            k_x0 = params.kappa(i, i, x, E) * x;
        }
        catch (std::invalid_argument &ex) {
            throw ex;
        }
        if (params.isOpen(i, E)) {
            em(i, i) = {cos(k_x0), -sin(k_x0)};
        } else {
            em(i, i) = cosh(k_x0);
        }
    }
    return em;
}

/*!
 * This method calculates \f$ \mathbf{E}^+ \f$ matrix for for a given point \f$x\f$ and given energy.
 * The matrix is diagonal and its elements are calculated the following way:
 * \f{itemize}{
 * \item $ \mathbf{E}^+_{n, n}(x, E) = \exp(i k x)$ if channel $n$ is open
 * \item $ \mathbf{E}^+_{n, n}(x, E) = \sinh(k x)$ if channel $n$ is closed
 * \item $ \mathbf{E}^+_{n, m}(x, E) = 0$ for $n\neq m$
 * \f}
 * @param x
 * @param E - energy
 * @return \f$ \mathbf{E}^+(x, E)\f$
 * @throws std::invalid_argument if x is wrong
 */
arma::cx_mat NonconstantGridSolver::calculateEP(double x, double E) {
    arma::cx_mat ep = params.Id();
    double k_x0;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        try {
            k_x0 = params.kappa(i, i, x, E) * x;
        }
        catch (std::invalid_argument &ex) {
            throw ex;
        }
        if (params.isOpen(i, E)) {
            ep(i, i) = {cos(k_x0), sin(k_x0)};
        } else {
            ep(i, i) = sinh(k_x0);
        }
    }
    return ep;
}

/*!
 * This method performs the Numerov iteration for a given energy for a case of some particular symmetry.
 *
 *  The initial value \f$ \mathbf{R}_0^{-1} \f$:
 *  \f{itemize}{
 *   \item $ \mathbf{R}_0^{-1} = \mathbf{0} $ if no symmetries
 *
 *   \item $ \mathbf{R}_0^{-1} = \mathbf{U}_0^{-1}(\mathbf{I} + \mathbf{B}) $ if the symmetry is described by
 *  $ \mathbf{B} $
 *  \f}
 *  \param[in] B - \f$ \mathbf{B}\f$ matrix to calculate the initial value
 *  \param[in] E - energy
 *  \return \f$ \mathbf{R}_{N-1}\f$
 *  \throws std::invalid_argument if \f$ \mathbf{U}_i\f$ cannot be calculated for given iteration \f$i\f$
 *  \throws std::runtime_error
 */

arma::cx_mat NonconstantGridSolver::fwdIteration(const arma::cx_mat &B, double E, std::vector<double> grid) {
    arma::cx_mat R0Inv, RBef, RBefBef;
    if (params.getNSymmetries() == 0) {
        R0Inv = arma::cx_mat(params.getNChannels(), params.getNChannels(), arma::fill::zeros);
    } else {
        try {
            R0Inv = arma::inv(U(grid[0], E, grid[1] - grid[0])) * (params.Id() + B); //R_0^{-1}
        } catch (std::invalid_argument &ex) {
            throw ex;
        } catch (...) {
            throw std::runtime_error("Error occured when calculating R0_INV for E = " + std::to_string(E));
        }
    }

    try {
        RBefBef = U(grid[1], E, grid[2] - grid[1]) - R0Inv;  //R_1
        RBef = U(grid[2], E, grid[3] - grid[2]) - arma::inv(RBefBef);
    } catch (std::invalid_argument &ex) {
        throw ex;
    } catch (...) {
        std::runtime_error("Error occured when calculating R1/R2 for E = " + std::to_string(E));
    }

    arma::cx_mat Ri, Ui;// R_{i}, U_{i}
    for (auto i = 3; i < grid.size() - 1; ++i) {
        try {
            auto dx_prev = grid[i] - grid[i - 1];
            auto dx = grid[i + 1] - grid[i];
            Ui = U(grid[i], E, dx);

            //if to double
            if (dx / dx_prev > 1.5)
                Ri = Ui - arma::inv(RBef * RBefBef);
                //if to halve
            else if (dx / dx_prev < 0.75) {
                arma::cx_mat T_alpha = T((grid[i] + grid[i - 1]) / 2, E, dx);
                Ri = Ui - 0.5 * (Id - T_alpha) * (arma::inv(Id + 2. * T_alpha) * (Id + arma::inv(RBef)) *
                                                  (Id - 4. * T(grid[i], E, dx))) * arma::inv(Id - T(grid[i], E, dx));
            }
                //if to leave the same
            else
                Ri = Ui - arma::inv(RBef);
        } catch (std::invalid_argument &ex) {
            throw ex;
        } catch (...) {
            throw std::runtime_error(
                    "Error occured when calculating R_" + std::to_string(i) + " for E = " + std::to_string(E));
        }
        RBefBef = RBef;
        RBef = Ri;
    }
    return Ri;
}
/*!
 * This method calculates the scattering matrix \f$ \mathbf{S}(E)\f$.
 * Its value is given by
 * \f{equation}{ \mathbf{S} = (\mathbf{R}_{N-1} \mathbf{e}^+_{N-1} - \mathbf{e}^+_{N}){-1} (\mathbf{R}_{N-1}\mathbf{e}^-_{N-1} - \mathbf{e}^-_{N}) \f}
 *
 * where \f$  \mathbf{e}^{\pm}_{i} =  (\mathbf{I} -  \mathbf{T}_{i}) \mathbf{E}^{\pm}_{i}\f$.
 * @param R_N - \f$ \mathbf{R}_{N-1}\f$
 * @param E - energy
 * @return \f$ \mathbf{S}\f$
 * @throws std::runtime_error if there is a problem with calculating
 */
arma::cx_mat NonconstantGridSolver::calculateS(const arma::cx_mat R_N, double E, std::vector<double> grid) {
    auto dx = grid[grid.size() - 1] - grid[grid.size() - 2];
    arma::cx_mat ITN, ITN1, epN, epN1, emN, emN1;

    ITN = params.Id() - T(grid[grid.size() - 2], E, dx);
    ITN1 = params.Id() - T(grid[grid.size() - 1], E, dx);
    epN = ITN * calculateEP(grid[grid.size() - 2], E);
    epN1 = ITN1 * calculateEP(grid[grid.size() - 1], E);
    emN = ITN * calculateEM(grid[grid.size() - 2], E);
    emN1 = ITN1 * calculateEM(grid[grid.size() - 1], E);
    //calculate S
    arma::cx_mat S;
    try {
        S = arma::inv(R_N * epN - epN1) * (R_N * emN - emN1);
    } catch (...) {
        throw std::runtime_error("Error when calculating S matrix for E = " + std::to_string(E));
    }

    return S;
}

void NonconstantGridSolver::saveS(const arma::cx_mat &S, const std::string S_type, const double E) {

    arma::mat reS = arma::real(S);
    arma::mat imS = arma::imag(S);

//    auto filenameReS = params.getDirID() + "S//reS" + S_type + std::to_string(E) + ".dat";
//    reS.save(filenameReS, arma::raw_ascii);
//
//    auto filenameImS = params.getDirID() + "S//imS" + S_type + std::to_string(E) + ".dat";
//    imS.save(filenameImS, arma::raw_ascii);

}

void NonconstantGridSolver::solveForEnergies(std::string directory) {
    //for(int i = 0; i < params.getNSymmetries(); i++)
    int i = 0;
    do {
        auto B = params.getB(i);
        std::vector<std::vector<double>> data(params.getNE());
        int j;
#pragma omp parallel for
//
        for (j = 0; j < params.getNE(); ++j) {
            auto E = params.getE(j);
            std::vector<double> row;
            try {
                std::vector<double> grid = generateGrid(E);
                arma::cx_mat R_N = fwdIteration(B, E, grid);

                arma::cx_mat S = calculateS(R_N, E, grid);
                auto S_value = S(0, 0);
                row.push_back(E);
                row.push_back(std::real(S_value));
                row.push_back(std::imag(S_value));
                data[j] = row;
////                saveS(S, "_sym=" + std::to_string(i) + "_E=", j, directory );
            }
            catch (std::runtime_error &ex) {
                std::cout << "Problem with solving for energy: " << E << '\n';
                throw ex;
            }
            catch (...) {
                std::cout << "ERROR in solveForEnergies\n";
            }
        }
        std::ofstream file(directory);
        for (auto row:data) {
            for (auto value:row) {
                file << value << " ";
            }
            file << "\n";
        }
        file.close();
        i++;
    } while (i < params.getNSymmetries());

}

/*!
 * This method generates grid points for Numerov calculations based on the energy E and potential.
 * @param E
 * @return x
 */
std::vector<double> NonconstantGridSolver::generateGrid(double E) {
    std::vector<double> x;
    double h = params.requiredDx(params.getXMin(), E);
    x = {params.getXMin(),
         params.getXMin() + h,
         params.getXMin() + 2 * h,
         params.getXMin() + 3 * h
    };
    auto dx_prev_prev = x[2] - x[1];
    auto dx_prev = x[3] - x[2];
    auto xValue = x[3];
    double dx;
    while (true) {
        auto req_dx = params.requiredDx(xValue, E);
        auto eps = 0.01;
        //halving
        if (req_dx < (1 - eps) * dx_prev) {
            dx = dx_prev / 2;
        }
        //doubling
        else if (req_dx > (1 + eps) * dx_prev and dx_prev == dx_prev_prev) {
            dx = 2 * dx_prev;
        }
        //the same
        else {
            dx = dx_prev;
        }
        xValue += dx;
        if (xValue > params.getXMax()) {
            break;
        }
        x.push_back(xValue);
        dx_prev_prev = dx_prev;
        dx_prev = dx;
    }
    return x;
}
