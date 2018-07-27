//
// Created by martas on 07.12.17.
//



#include "ConstantGridSolver.h"
/*!
 * \file
 * \brief Definitions of Solver class methods.
 *
 */

/*!
 * This method calculates U matrix for given index j and set of parameters
 * \param[in] j - index of x value for which U is calculated
 * \param[in] params - Parameters object for which U is calculated
 * \return U matrix
*/
arma::cx_mat ConstantGridSolver::calculateU(int j, double E) {
    arma::cx_mat T,  IT;
    try{
        T= calculateT(j, E);
    }
    catch(std::invalid_argument &ex){
        throw ex;
    }
      IT = params.Id() - T;
    return 12. * arma::inv(IT) - 10. * params.Id();
}

/*!
 * This method calculates and stores the sets of matrices R (as \link Solver::RList RList\endlink), F (\link Solver::FList
 * PsiList\endlink), Psi (\link Solver::PsiList PsiList\endlink) and sets of vectors f and psi (\link Solver::psiList psiList\endlink).
 * The results are stored (as the attributes of the Solver object), but are not yet saved to file.
 * \param[in] &params - Parameters object to perform calculations for.
 */

arma::cx_mat ConstantGridSolver::fwdIteration(const arma::cx_mat &B, double E) {

    arma::cx_mat R0Inv, RBef;
    try {
        R0Inv = arma::inv(calculateU(0, E)) * (params.Id() + B); //R_0^{-1}
    }catch(std::invalid_argument &ex){
        throw ex;
    }catch(...){
        throw std::runtime_error("Error occured when calculating R0_INV for E = " + std::to_string(E));
    }
    try{
        RBef = calculateU(1, E) - R0Inv;  //R_1
    }catch(std::invalid_argument &ex){
        throw ex;
    }catch(...){
        std::runtime_error("Error occured when calculating R1 for E = " + std::to_string(E));
    }

    arma::cx_mat Ri, Ui;// R_{i}, U_{i}
    for (auto i = 2; i < params.NX() - 1; ++i) {
        try {
            Ui = calculateU(i, E);
            Ri = Ui - arma::inv(RBef);
        }catch(std::invalid_argument &ex){
            throw ex;
        }catch(...){
            throw std::runtime_error("Error occured when calculating R_" +std::to_string(i) + " for E = " + std::to_string(E));
        }
        RBef = Ri;
    }
    return Ri;
}

arma::cx_mat ConstantGridSolver::calculateS(const arma::cx_mat R_N, double E) {

    arma::cx_mat ITN = params.Id() - calculateT(-2, E);
    arma::cx_mat ITN1 = params.Id() - calculateT(-1, E);

    arma::cx_mat epN = ITN * calculateEP(-2, E);
    arma::cx_mat epN1 = ITN1 * calculateEP(-1, E);
    arma::cx_mat emN =  ITN * calculateEM(-2, E);
    arma::cx_mat emN1 = ITN1 * calculateEM(-1, E);

    //calculate S
    arma::cx_mat S;
    try{
        S = arma::inv(R_N * epN - epN1) * (R_N * emN - emN1);
    }catch(...){
        throw std::runtime_error("Error when calculating S matrix for E = " + std::to_string(E));
    }

    return S;
}

void ConstantGridSolver::saveS(const arma::cx_mat &S, const std::string S_type, const double E, const std::string directory) {

    arma::mat reS = arma::real(S);
    arma::mat imS = arma::imag(S);

    auto filenameReS = directory + "re_S" + S_type + std::to_string(E) + ".dat";
    reS.save(filenameReS, arma::raw_ascii);

    auto filenameImS = directory + "im_S" + S_type + std::to_string(E) + ".dat";
    imS.save(filenameImS, arma::raw_ascii);

}

/*!
 * This method calculates T matrix for given index j and set of parameters
 * \params[in] j - index of x value for which T is calculated
 * \params[in] params - Parameters object for which T is calculated
 * \return T matrix
 */
arma::cx_mat ConstantGridSolver::calculateT(int j, double E) const{
    arma::cx_mat V;
    try{
        V = params.getVMatrix(j);
    }
    catch(std::invalid_argument &ex){
        throw ex;
    }
    return -(1. / params.getUnit()) * params.getDx() * params.getDx() * (E * params.Id() - V) / 12.;
}


void ConstantGridSolver::modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E) {
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        if (!params.isOpen(i, E)) {
            j1(i, i) = j1(i, i) / j0(i, i);
            j0(i, i) = 1;
            n1(i, i) = n1(i, i) / n0(i, i);
            n0(i, i) = 1;
        }
    }
}

ConstantGridSolver::ConstantGridSolver(const Parameters& params) : params(params) {
    ;
    std::cout << "ConstantGridSolver done\n";
}

arma::cx_mat ConstantGridSolver::calculateEM(int ind_x, double E) {
    arma::cx_mat em = params.Id();
    double k_x0;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        try{
            k_x0 = params.kappa(i, i, ind_x, E) * params.x(ind_x);
        }
        catch(std::invalid_argument &ex){
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

arma::cx_mat ConstantGridSolver::calculateEP(int ind_x, double E) {
    arma::cx_mat ep = params.Id();
    double k_x0;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        try{
            k_x0 = params.kappa(i, i, ind_x, E) * params.x(ind_x);
        }
        catch(std::invalid_argument &ex){
            throw ex;
        }
        if (params.isOpen(i, E)) {
            ep(i, i) = {cos(k_x0), sin(k_x0)};
        } else{
            ep(i, i) =sinh(k_x0);
        }
    }
    return ep;
}

void ConstantGridSolver::solveForEnergies(std::string directory) {

    for(int i = 0; i < params.getNSymmetries(); i++){
        auto B = params.getB(i);
        for (int j = 0; j < params.getNE(); ++j) {
            auto E = params.getE(j);
            try{
                auto R_N = fwdIteration(B, E);
                auto S = calculateS(R_N, E);
                saveS(S, "_sym=" + std::to_string(i) + "_E=", E, directory );
            }
            catch(std::runtime_error &ex){
                std::cout << "Problem with solving for energy: " << E << '\n';
                throw ex;
            }
        }
    }

}


































