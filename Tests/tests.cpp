//
// Created by martas on 22.07.18.
//

#include <gtest/gtest.h>
#include "ParametersTests.h"
#include "ConstantGridSolverTests.h"
#include "NonconstantGridSolverTests.h"


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}