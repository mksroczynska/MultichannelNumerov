//
// Created by martas on 09.08.18.
//

#ifndef NUMEROV_NONCONSTANTGRIDSOLVERTESTS_H
#define NUMEROV_NONCONSTANTGRIDSOLVERTESTS_H

#include <gtest/gtest.h>
#include "../NonconstantGridSolver.h"
#include "ParametersTests.h"

class NonconstantGridSolver_Test: public ::testing::Test{
public:
    NonconstantGridSolver solver;

    void SetUp() override {
        const std::vector<std::string> filenamesList= {
                testing_data_dir + "good_file_with_params_nx=4.txt",
                testing_data_dir + "good_EList_doubles.dat",
                testing_data_dir + "good_V_two_channels_nx=4_to_nonconstgrid.txt",
                testing_data_dir + "good_B_two_channels"
        };
        const Parameters parameters_testing(filenamesList);
        solver.setParameters(parameters_testing);
    }
};
class NonconstantGridSolver_Q_Test: public NonconstantGridSolver_Test{};
TEST_F(NonconstantGridSolver_Q_Test, returnsCorrectValueForCorrectParameters){
    arma::mat expected_re = {{4, -2}, {-3, 1}};
    arma::mat expected_im(2,2,arma::fill::zeros);

    arma::cx_mat Q = solver.Q(0., 5);

    EXPECT_TRUE(arma::approx_equal(expected_re, arma::real(Q), "absdiff", 0.0000001));
    EXPECT_TRUE(arma::approx_equal(expected_im, arma::imag(Q), "absdiff", 0.0000001));
}
class NonconstantGridSolver_generateGrid_Test: public NonconstantGridSolver_Test{};
TEST_F(NonconstantGridSolver_generateGrid_Test, returnsCorrrectGridForCorrectValues){
    auto grid = solver.generateGrid(15);
    for(auto x: grid){
        std::cout << x << ", ";
    }
}

#endif //NUMEROV_NONCONSTANTGRIDSOLVERTESTS_H
