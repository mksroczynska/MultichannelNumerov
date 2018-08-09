//
// Created by martas on 22.06.18.
//

#include "NonconstantGridSolver.h"

arma::cx_mat NonconstantGridSolver::calculateT(double x, double E, double dx) {
    return -(1. / params.getUnit()) * dx * dx * (E * Id - params.getV(x)) / 12.;
}

arma::cx_mat NonconstantGridSolver::calculateU(double x, double E, double dx) {
    auto T = calculateT(x, E, dx);//calculate T
    return 12. * arma::inv(Id - T) - 10. * Id;
}

arma::cx_mat
NonconstantGridSolver::calculateFAlpha(double x, double dx, double E, const arma::cx_mat &F_i, const arma::cx_mat &F_i1,
                                    const arma::cx_mat &U_i) {
    auto T_alpha = calculateT(x, E, dx);

    return (Id - T_alpha) * arma::inv(2. * Id + T_alpha) * ((Id + U_i) * F_i - F_i1);

}

void NonconstantGridSolver::modifyCCnj(arma::cx_mat &n1, arma::cx_mat &n0, arma::cx_mat &j1, arma::cx_mat &j0, double E) {
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        if (!params.isOpen(i, E)) {
            j1(i, i) = j1(i, i) / j0(i, i);
            j0(i, i) = 1;
            n1(i, i) = n1(i, i) / n0(i, i);
            n0(i, i) = 1;
        }
    }
}
auto NonconstantGridSolver::calculateEMP(double x0, double x1, double E) {
    arma::cx_mat Zeros(n, n, arma::fill::zeros);
    auto em0 = Zeros;
    auto em1 = Zeros;
    auto ep0 = Zeros;
    auto ep1 = Zeros;
    double k_x0, k_x1;
    for (unsigned int i = 0; i < params.getNChannels(); ++i) {
        k_x0 = params.kappa(i, i, x0, E) * x0;
        k_x1 = params.kappa(i, i, x1, E) * x1;
        if (params.isOpen(i, E)) {
            em0(i, i) = {cos(k_x0), -sin(k_x0)};
            em1(i, i) = {cos(k_x1), -sin(k_x1)};
            ep0(i, i) = {cos(k_x0), sin(k_x0)};
            ep1(i, i) = {cos(k_x1), sin(k_x1)};
        }
        else {
            //em0(i, i) = exp(-k_x0);
            //em1(i, i) = exp(-k_x1);
            em0(i, i) = cosh(k_x0);
            em1(i, i) = cosh(k_x1);
            ep0(i, i) = sinh(k_x0);
            ep1(i, i) = sinh(k_x1);
        }
    }
    return std::make_tuple(em0, em1, ep0, ep1);
}

/*!
 * This method calculates and stores the sets of matrices R (as \link Solver::RList RList\endlink), F (\link Solver::FList
 * PsiList\endlink), Psi (\link Solver::PsiList PsiList\endlink) and sets of vectors f and psi (\link Solver::psiList psiList\endlink).
 * The results are stored (as the attributes of the Solver object), but are not yet saved to file.
 * \param[in] &params - Parameters object to perform calculations for.
 */

arma::cx_mat NonconstantGridSolver::fwdIteration(const arma::cx_mat &B, double E) {
    //initial grid: zalezny od potencjalu w xmin
    auto dxBef = params.requiredDx(params.getXMin(), E);
    std::cout << "dxBef: " << dxBef <<'\n';
    arma::cx_mat R0Inv = arma::inv(calculateU(params.getXMin(), E, dxBef)) * (Id + B); //R_0^{-1}
    std::cout << "R0Inv: ";
    R0Inv.print();
    arma::cx_mat RBefBef = calculateU(params.getXMin()+dxBef, E, dxBef) - R0Inv; //R_1
    std::cout << "R_1: ";
    RBefBef.print();
    auto xBef = params.getXMin() + 2 * dxBef;
    arma::cx_mat RBef = calculateU(xBef, E, dxBef) - arma::inv(RBefBef);  // initial value of R_2
    std::cout << "R_2: ";
    RBef.print();
    auto dx = changeDxIfNeeded(xBef, dxBef, E);
    std::cout << "New dx: " << dx << '\n';
//
   arma::cx_mat Ri, Ui, Fi;// R_{i}, U_{i}
//
    if(gridPoints.empty()){
        gridPoints.push_back(params.getXMin());
        gridPoints.push_back(params.getXMin()+dxBef);
        gridPoints.push_back(xBef);
    }

    // calculate R up to N (R_N)
    for (auto x = xBef + dx ; x < params.getXMax(); x += dx) {
////        std::cout << "FwdIteration: x=" <<x << '\n';
//            gridPoints.push_back(x);
//            Ui = calculateU(x, E, dx);
//
//        //sprawdzac jakos czy zmieniony i jak
//
//        //the same
//        if(dxBef == dx) {
//            Ri = Ui - arma::inv(RBef); // calculate R_{i}
//            Fi = RBef * FBef;          // F_{i} = R_{i-1} * F_{i-1}
//        }
//
//        //double grid
//        else if(dxBef < dx){
//            Ri = Ui - arma::inv(RBef * RBefBef); // R_{i} w punkcie oddalonym o wieksze dx od poprzedniego
//            Fi = Ri * FBefBef;  // R_{i} = F_{i}/F_{i-2} (i-2 jest w punkcie oddalonym o powiekszone dx od x_{i})
//
//        }
//
//        //halving
//        else if(dxBef > dx){
//            Fi = Ui*FBef - FBefBef; // F_{i} w punkcie oddalonym o wcześniejszy grid od poprzedniego punktu}
//            arma::cx_mat Falpha = calculateFAlpha(x, dx, E, Fi, FBef, Ui);
//            Ri = Ui - Falpha * arma::inv(FBef) * arma::inv(RBef);
//            Fi = Falpha;
////        }
//
//
//        //calculate F
//        //gdzie obliczac dx?
//        //trzymac wszystkie x w vectorze?
//        //trzymac wszystkie F, R?
//        //zrobic vector vectorow trzymajacy wszystkie rozwiazania?
//        // + metody zapisujace do pliku
//        //przetestowac na 1x1
//
//        RBefBef = RBef;
//        RBef = Ri;
//            FBefBef = FBef;
//            FBef = Fi;
//        dxBef = dx;
//        dx = changeDxIfNeeded(x, dxBef, E);
    }
    return Ri; // R_{N} = Ri
}

arma::cx_mat NonconstantGridSolver::calculateS(const arma::cx_mat R_N, double E, double dx) {
    //define EM_{N} and EM_{N+1}
    auto EMP = calculateEMP(params.getXMax(), params.getXMax()+ dx, E);

    // calculate ep_{N}, ep_{N+1} and em_{N}, em_{N+1}
    auto TN = calculateT(params.getXMax(), E, dx);
    auto TN1 = calculateT(params.getXMax() + dx, E, dx);

    auto epN = (Id - TN) * std::get<2>(EMP);
    auto epN1 = (Id - TN1) * std::get<3>(EMP);
    auto emN = (Id - TN) * std::get<0>(EMP);
    auto emN1 = (Id - TN1) * std::get<1>(EMP);

    //calculate S
    return arma::inv(R_N * epN - epN1) * (R_N * emN - emN1);
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

double NonconstantGridSolver::changeDxIfNeeded(double xBef, double dxBef, double E) {

    auto dx_req = params.requiredDx(xBef, E);

    if (abs(dx_req - dxBef) <0.1*abs(dxBef)) return dxBef; //the same grid
    else if (dx_req < dxBef) return dxBef / 2.; //halving
    else return 2. * dxBef; //doubling
}

void NonconstantGridSolver::solveForEnergies() {

    std::cout << "Solving ...";

    auto B0 = params.getB(0);
    auto B1 = params.getB(1);
    fwdIteration(B1, params.getE(35));

//    int j;
//    for (j = 0; j < params.getNE(); ++j) {
//
//        auto E = params.getE(j);
//
//        auto R_N_p = fwdIteration(B0, E);
////        auto S_p = calculateS(R_N_p, E);
////        saveS(S_p, "p", E);
//
//        auto R_N_m = fwdIteration(B1, E);
////        auto S_m = calculateS(R_N_m, E);
////        saveS(S_m, "m", E);
//
//    }
    std::cout << "Final number of gridpoints: " << gridPoints.size() << '\n';
    std::ofstream plik("xPoints.dat");
    for (auto gridPoint : gridPoints) {
        plik << gridPoint;
        plik << "\n";
    }
    plik.close();
}