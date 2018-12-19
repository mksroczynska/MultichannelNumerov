//
// Created by martas on 07.12.17.
//



#include "ConstantGridSolver.h"
/*!
 * \file
 * \brief Definitions of ConstantGridSolver class methods.
 *
 */

/*!
 * This method calculates \f$ \mathbf{U} \f$ matrix for given index j using the set of parameters provided to the ConstantGridSolver
 * object according to the following formula:
 * \f{equation}{
 * \mathbf{U}_j = 12 (\mathbf{I} - \mathbf{T}_j)^{-1} - 10 \mathbf{I}. \f}
 * \param[in] j - index on the grid of x value
 * \param[in] E - energy value
 * \return \f$ \mathbf{U}(x_j, E) \f$
 * \throws std::invalid_argument if \f$ x_j\f$ does not exist
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
 * This method performs the Numerov iteration for a given energy for a case of some particular symmetry.
 *
 *  The initial value \f$ \mathbf{R}_0^{-1} \f$:
 *  \f{itemize}{
 *   \item $ \mathbf{R}_0^{-1} = \mathbf{0} $ if no symmetries
 *
 *   \item $ \mathbf{R}_0^{-1} = \mathbf{U}_0^{-1}(\mathbf{I} + \mathbf{B}) $ if the symmetry is described by
 *  $ \mathbf{B} $
 *  \f}
 *
 *  Every value depends on the previous one: \f$ \mathbf{R}_j = \mathbf{U}_j - \mathbf{R}_{j-1}^{-1}. \f$
 *  \param[in] B - \f$ \mathbf{B}\f$ matrix to calculate the initial value
 *  \param[in] E - energy
 *  \return \f$ \mathbf{R}_{N-1}\f$
 *  \throws std::invalid_argument if \f$ \mathbf{U}_i\f$ cannot be calculated for given iteration \f$i\f$
 *  \throws std::runtime_error
 */

arma::cx_mat ConstantGridSolver::fwdIteration(const arma::cx_mat &B, double E) {

    arma::cx_mat R0Inv, RBef;
    if(params.getNSymmetries() == 0){
        R0Inv  = arma::cx_mat(params.getNChannels(), params.getNChannels(), arma::fill::zeros);
    }
    else{
        try {
            R0Inv = arma::inv(calculateU(0, E)) * (params.Id() + B); //R_0^{-1}
        }catch(std::invalid_argument &ex){
            throw ex;
        }catch(...){
            throw std::runtime_error("Error occured when calculating R0_INV for E = " + std::to_string(E));
        }
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
/*!
 * This method saves the scattering matrix in a given directory.
 * The real and imaginary part of \f$ \mathbf{S} \f$ are saved in separate files.
 *
 * Paths:
 *
 * \f$ Re(\mathbf{S}) \f$: directory/re_SE.dat (E is the value of the energy)
 * \f$ Im(\mathbf{S}) \f$: directory/im_SE.dat (E is the value of the energy)
 *
 * @param S - scattering matrix to be saved
 * @param E - energy
 * @param directory - where to save the files
 */
void ConstantGridSolver::saveS(const arma::cx_mat &S, const int E, const std::string directory) {

    arma::mat reS = arma::real(S);
    arma::mat imS = arma::imag(S);

    auto filenameReS = directory + "re_S" + std::to_string(E) + ".dat";
    reS.save(filenameReS, arma::raw_ascii);

    auto filenameImS = directory + "im_S" + std::to_string(E) + ".dat";
    imS.save(filenameImS, arma::raw_ascii);

}

/*!
 * This method calculates \f$ \mathbf{T}\f$ matrix for a given point \f$x_j\f$ on the grid and given energy
 * according to the formula:
 * \f{equation}{
 * \mathbf{T}_j = -\frac{dx}{12}\mathbf{Q}_j.
 *  \f}
 * @param j - index of the value on the grid
 * @param E - energy
 * @return \f$ \mathbf{T}(x_j, E)\f$
 * @throws std::invalid_argument if the index \f$ j \f$ is wrong
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
}

/*!
 * This method calculates \f$ \mathbf{E}^- \f$ matrix for for a given point \f$x_j\f$ on the grid and given energy.
 * The matrix is diagonal and its elements are calculated the following way:
 * \f{itemize}{
 * \item $ \mathbf{E}^-_{n, n}(x_j, E) = \exp(-i k x_j)$ if channel $n$ is open
 * \item $ \mathbf{E}^-_{n, n}(x_j, E) = \cosh(k x_j)$ if channel $n$ is closed
 * \item $ \mathbf{E}^-_{n, m}(x_j, E) = 0$ for $n\neq m$
 * \f}
 * @param j - index of the value on the grid
 * @param E - energy
 * @return \f$ \mathbf{E}^-(x_j, E)\f$
 * @throws std::invalid_argument if j is wrong
 */
arma::cx_mat ConstantGridSolver::calculateEM(int j, double E) {
    arma::cx_mat em = params.Id();
    double k_x0;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        try{
            k_x0 = params.kappa(i, i, j, E) * params.x(j);
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

/*!
 * This method calculates \f$ \mathbf{E}^+ \f$ matrix for for a given point \f$x_j\f$ on the grid and given energy.
 * The matrix is diagonal and its elements are calculated the following way:
 * \f{itemize}{
 * \item $ \mathbf{E}^+_{n, n}(x_j, E) = \exp(i k x_j)$ if channel $n$ is open
 * \item $ \mathbf{E}^+_{n, n}(x_j, E) = \sinh(k x_j)$ if channel $n$ is closed
 * \item $ \mathbf{E}^+_{n, m}(x_j, E) = 0$ for $n\neq m$
 * \f}
 * @param j - index of the value on the grid
 * @param E - energy
 * @return \f$ \mathbf{E}^+(x_j, E)\f$
 * @throws std::invalid_argument if j is wrong
 */
arma::cx_mat ConstantGridSolver::calculateEP(int j, double E) {
    arma::cx_mat ep = params.Id();
    double k_x0;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        try{
            k_x0 = params.kappa(i, i, j, E) * params.x(j);
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

    //for(int i = 0; i < params.getNSymmetries(); i++)
    int i = 0;
    do {
        auto B = params.getB(i);
        std::vector<std::string> data(params.getNE());
        int j;
        #pragma omp parallel for

        for ( j = 0; j < params.getNE(); ++j) {
            auto E = params.getE(j);
            try{
                auto R_N = fwdIteration(B, E);
                auto S = calculateS(R_N, E);
                data[j] = tostring_E_S(S, E);
                //saveS(S, "_sym=" + std::to_string(i) + "_E=", j, directory );
            }
            catch(std::runtime_error &ex){
                std::cout << "Problem with solving for energy: " << E << '\n';
                throw ex;
            }
        }
        std::ofstream file(directory);
        for(auto row:data){
            file << row;
            file << "\n";
        }
        file.close();
        i++;
    }while(i < params.getNSymmetries());

}

/*!
 * Single row for the file with results.
 * Format: E re(S(0,0)) re(S(0,1)) ... re(S(n, 0)... re(S(n, n)) im(S(0,0)) im(S(0,1)) ... im(S(n, 0)... im(S(n, n))
 * @param S
 * @param E
 * @return
 */
std::string ConstantGridSolver::tostring_E_S(arma::cx_mat S, double E) {
    std::string S_re = "";
    std::string S_im = "";
    int n = std::sqrt(S.size());
    for (int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j){
           S_re += std::to_string(std::real(S(i, j)));
           S_re += " ";
           S_im += std::to_string(std::imag(S(i, j)));
           S_im += " ";
        }
    }
    return std::to_string(E) + " " + S_re + S_im;
}


































