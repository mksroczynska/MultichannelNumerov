#ifndef NUMEROV_CONSTANGRIDSOLVERTESTS_H
#define NUMEROV_CONSTANGRIDSOLVERTESTS_H

#include <gtest/gtest.h>
#include "../ConstantGridSolver.h"


class ConstantGridSolver_Test: public ::testing::Test{
public:
    ConstantGridSolver solver;
    
    void SetUp() override {
        const std::vector<std::string> filenamesList= {
                testing_data_dir + "good_file_with_params_nx=4.txt",
                testing_data_dir + "good_EList_doubles.dat",
                testing_data_dir + "good_V_two_channels_nx=4.txt",
                testing_data_dir + "good_B_two_channels"
        };
        const Parameters parameters_testing(filenamesList);
        solver.setParameters(parameters_testing);
    }
};
class ConstantGridSolver_calculateT_Test: public ConstantGridSolver_Test{};
class ConstantGridSolver_calculateU_Test: public ConstantGridSolver_Test{};
class ConstantGridSolver_calculateEP_Test: public ConstantGridSolver_Test{};
class ConstantGridSolver_calculateEM_Test: public ConstantGridSolver_Test{};
class ConstantGridSolver_fwdIteration_Test: public ConstantGridSolver_Test{};
class ConstantGridSolver_calculateS_Test: public ConstantGridSolver_Test{};
class ConstantGridSolver_solveForEnergies_Test: public ConstantGridSolver_Test{};
//******************calculateT(i, E, Id)****************************

TEST_F(ConstantGridSolver_calculateT_Test, worksCorrectForCorrectArguments){
    arma::cx_mat T = solver.calculateT(0, 5);
    arma::mat expected_re = {{-(0.25/12) * 4., (0.25/12 )* 2},
                                {0.25/12 * 3, -0.25/12}};
    arma::mat expected_im = arma::zeros(2,2);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(T(i, j)));
            EXPECT_DOUBLE_EQ(expected_im(i, j), std::imag(T(i,j)));
        }
    }

}
TEST_F(ConstantGridSolver_calculateT_Test, worksCorrectForCorrectArgumentsWithNegativeIndex){
    arma::cx_mat T = solver.calculateT(-4, 5);
    arma::mat expected_re = {{-(0.25/12) * 4., (0.25/12 )* 2},
                             {0.25/12 * 3, -0.25/12}};
    arma::mat expected_im = arma::zeros(2,2);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(T(i, j)));
            EXPECT_DOUBLE_EQ(expected_im(i, j), std::imag(T(i,j)));
        }
    }

}
TEST_F(ConstantGridSolver_calculateT_Test, throwsInvalidArgumentExceptionForIncorrectIndexJ){
    EXPECT_THROW(solver.calculateT(10, 5), std::invalid_argument);
    EXPECT_THROW(solver.calculateT(-10, 5), std::invalid_argument);
}



//************************U(i, E)*******************

TEST_F(ConstantGridSolver_calculateU_Test, worksCorrectForCorrectArguments){
    arma::cx_mat U = solver.calculateU(0, 5);
    arma::mat expected_re =
            {{1.1030684500393395, 0.4531864673485445},
            {0.6797797010228168, 1.7828481510621579}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(U(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(U(i,j)));
        }
    }
}

TEST_F(ConstantGridSolver_calculateU_Test, rethrowsExceptionFromCalculateTIfCatches){
    EXPECT_THROW(solver.calculateU(-10, 5), std::invalid_argument);
}

//************************calculateEP(i, E)******************
TEST_F(ConstantGridSolver_calculateEP_Test, worksCorrectForCorrectValuesTwoOpenChannels){

    EXPECT_NO_THROW(solver.calculateEP(0, 13));
    arma::cx_mat EP = solver.calculateEP(0,13);
    arma::mat expected_re =
            {{1, 0},
            {0, 1}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EP(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EP(i,j)));
        }
    }
    EXPECT_NO_THROW(solver.calculateEP(1, 13));
    arma::cx_mat EP2 = solver.calculateEP(1,13);
    expected_re =
            {{std::cos(std::sqrt(2.)), 0},
             {0, std::cos(0.5*std::sqrt(5.))}};
    arma::mat expected_im = {{std::sin(std::sqrt(2.)), 0}, {0, std::sin(0.5*std::sqrt(5.))}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EP2(i, j)));
            EXPECT_DOUBLE_EQ(expected_im(i,j), std::imag(EP2(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateEP_Test, worksCorrectForCorrectValuesTwoClosedChannels){

    EXPECT_NO_THROW(solver.calculateEP(0, 5));
    arma::cx_mat EP = solver.calculateEP(0,5);
    arma::mat expected_re =
            {{std::sinh(0), 0},
             {0, std::sinh(0)}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EP(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EP(i,j)));
        }
    }
    EXPECT_NO_THROW(solver.calculateEP(1, 5));
    arma::cx_mat EP2 = solver.calculateEP(1,5);
    expected_re =
            {{std::sinh(0), 0},
             {0, std::sinh(0.5*std::sqrt(3.))}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EP2(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EP2(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateEP_Test, worksCorrectForCorrectValuesOneOpenOneClosedChannel){
    EXPECT_NO_THROW(solver.calculateEP(0, 10));
    arma::cx_mat EP = solver.calculateEP(0,10);
    arma::mat expected_re =
            {{1, 0},
             {0, std::sinh(0)}};

    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EP(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EP(i,j)));
        }
    }
    EXPECT_NO_THROW(solver.calculateEP(1, 10));
    arma::cx_mat EP2 = solver.calculateEP(1,10);
    expected_re =
            {{std::cos(0.5 * std::sqrt(5.)), 0},
             {0, std::sinh(1./std::sqrt(2.))}};
    arma::mat expected_im = {{std::sin(0.5 * std::sqrt(5.)), 0}, {0, 0}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EP2(i, j)));
            EXPECT_DOUBLE_EQ(expected_im(i,j), std::imag(EP2(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateEP_Test, rethrowsExceptionFromKappaIfCatches){
    EXPECT_THROW(solver.calculateEP(-10, 5), std::invalid_argument);
}

//************************calculateEM(i, E)******************
TEST_F(ConstantGridSolver_calculateEM_Test, worksCorrectForCorrectValuesTwoOpenChannels){

    EXPECT_NO_THROW(solver.calculateEM(0, 13));
    arma::cx_mat EM = solver.calculateEM(0,13);
    arma::mat expected_re =
            {{1, 0},
             {0, 1}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EM(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EM(i,j)));
        }
    }
    EXPECT_NO_THROW(solver.calculateEM(1, 13));
    arma::cx_mat EM2 = solver.calculateEM(1,13);
    expected_re =
            {{std::cos(std::sqrt(2)), 0},
             {0, std::cos(0.5*std::sqrt(5.))}};
    arma::mat expected_im = {{-std::sin(std::sqrt(2)), 0}, {0, -std::sin(0.5*std::sqrt(5.))}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EM2(i, j)));
            EXPECT_DOUBLE_EQ(expected_im(i,j), std::imag(EM2(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateEM_Test, worksCorrectForCorrectValuesTwoClosedChannels){

    EXPECT_NO_THROW(solver.calculateEM(0, 5));
    arma::cx_mat EM = solver.calculateEM(0,5);
    arma::mat expected_re =
            {{std::cosh(0), 0},
             {0, std::cosh(0)}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EM(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EM(i,j)));
        }
    }
    EXPECT_NO_THROW(solver.calculateEM(1, 5));
    arma::cx_mat EM2 = solver.calculateEM(1,5);
    expected_re =
            {{std::cosh(0), 0},
             {0, std::cosh(0.5*std::sqrt(3.))}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EM2(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EM2(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateEM_Test, worksCorrectForCorrectValuesOneOpenOneClosedChannel){
    EXPECT_NO_THROW(solver.calculateEM(0, 10));
    arma::cx_mat EM = solver.calculateEM(0,10);
    arma::mat expected_re =
            {{1, 0},
             {0, std::cosh(0)}};

    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EM(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(EM(i,j)));
        }
    }
    EXPECT_NO_THROW(solver.calculateEM(1, 10));
    arma::cx_mat EM2 = solver.calculateEM(1,10);
    expected_re =
            {{std::cos(0.5 * std::sqrt(5.)), 0},
             {0, std::cosh(1/std::sqrt(2.))}};
    arma::mat expected_im = {{-std::sin(0.5 * std::sqrt(5.)), 0}, {0, 0}};
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_DOUBLE_EQ(expected_re(i,j), std::real(EM2(i, j)));
            EXPECT_DOUBLE_EQ(expected_im(i,j), std::imag(EM2(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateEM_Test, rethrowsExceptionFromKappaIfCatches){
    EXPECT_THROW(solver.calculateEM(-10, 5), std::invalid_argument);
}

//************************calculateS(R_N, E)******************
TEST_F(ConstantGridSolver_calculateS_Test, worksCorrectForCorrectValues){
    arma::mat R = {{4.239156504661 , 3.023409142108},
                    {3.2691144155848, 4.9762723250913}};
    arma::cx_mat RN(R, arma::zeros(2,2));
    arma::mat expected_re =
            {{1.012148643689, 0.005671732038210},
             {0.003902015446592, 1.0011947023937}};
    arma::cx_mat s;
    EXPECT_NO_THROW(s = solver.calculateS(RN, 5));
    for(int i = 0; i < 2; i++){
        for(int j = 0; j<2; j++){
            EXPECT_TRUE(std::abs(expected_re(i,j)- std::real(s(i, j))< 0.000001));
            EXPECT_DOUBLE_EQ(0, std::imag(s(i,j)));
        }
    }
}
TEST_F(ConstantGridSolver_calculateS_Test, throwsExceptionForIncorrectRN){
    arma::cx_mat RN(1,1, arma::fill::zeros);
    EXPECT_THROW(solver.calculateS(RN, 5), std::runtime_error);
}


//************************saveS(i, E)******************
//TEST_F(ConstantGridSolver_saveS_Test, savesUnderTheCorrectNames);
//TEST_F(ConstantGridSolver_saveS_Test, throwsExceptionIfCannotSave);

//************************fwdIteration(B, E)********************
TEST_F(ConstantGridSolver_fwdIteration_Test, returnsCorrectValueForCorrectArguments) {
    arma::mat expected_re =
            {{4.239156504661011,  3.023409142107999},
             {3.2691144155847773, 4.9762723250913465}};
    arma::cx_mat B(2, 2, arma::fill::eye);
    arma::cx_mat result;
    EXPECT_NO_THROW(result = solver.fwdIteration(B, 5));
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            EXPECT_DOUBLE_EQ(expected_re(i, j), std::real(result(i, j)));
            EXPECT_DOUBLE_EQ(0, std::imag(result(i, j)));
        }
    }
}

TEST_F(ConstantGridSolver_fwdIteration_Test, throwsExceptionIfCannotComputeR0Inv){
    arma::cx_mat B(3, 3, arma::fill::eye);
    EXPECT_THROW(solver.fwdIteration(B, 5), std::runtime_error);
}

//*************************solveForEnergies()*******************
TEST_F(ConstantGridSolver_solveForEnergies_Test, doesNotThrowIfEverythingGood){
solver.solveForEnergies("//home//martas//Dropbox//Doktoratowe//AIcollisions//MultichannelNumerovAIWithTests_1//MultichannelNumerovAIWithTests//Tests//TestData//");

    EXPECT_NO_THROW(solver.solveForEnergies("//home//martas//Dropbox//Doktoratowe//AIcollisions//MultichannelNumerovAIWithTests_1//MultichannelNumerovAIWithTests//Tests//TestData//"));
}

#endif