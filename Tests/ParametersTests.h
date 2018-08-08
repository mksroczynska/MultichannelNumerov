//
// Created by martas on 22.07.18.
//

#ifndef NUMEROV_PARAMETERSTESTS_H
#define NUMEROV_PARAMETERSTESTS_H

#include <gtest/gtest.h>
#include "../Parameters.h"

std::string testing_data_dir = "//home//martas//Dropbox//Doktoratowe//AIcollisions//MultichannelNumerovAIWithTests_1//MultichannelNumerovAIWithTests//Tests//TestData//";

class Parameters_loadParams_Test: public ::testing::Test {
    public:
    Parameters parameters;
};
class Parameters_loadV_Test: public ::testing::Test {
public:
Parameters parameters;
};

//******loadParams()*********

TEST_F(Parameters_loadParams_Test, worksGoodForGoodFile){
    EXPECT_NO_THROW(parameters.loadParams(testing_data_dir + "good_file_with_params.txt"));
    EXPECT_DOUBLE_EQ(0, parameters.getXMin());
    EXPECT_DOUBLE_EQ(1, parameters.getXMax());
    EXPECT_DOUBLE_EQ(0.5, parameters.getDx());
    EXPECT_EQ(1, parameters.getUnit());
    EXPECT_EQ(1, parameters.getNSymmetries());
    EXPECT_EQ(20, parameters.getGrid_points_per_lambda());

}
TEST_F(Parameters_loadParams_Test, worksGoodForGoodFileWithSomeDoubles){
    EXPECT_NO_THROW(parameters.loadParams(testing_data_dir + "good_file_with_params_some_doubles.txt"));
    EXPECT_DOUBLE_EQ(0, parameters.getXMin());
    EXPECT_DOUBLE_EQ(10, parameters.getXMax());
    EXPECT_DOUBLE_EQ(0.1, parameters.getDx());
    EXPECT_EQ(1, parameters.getUnit());
    EXPECT_EQ(1, parameters.getNSymmetries());
    EXPECT_EQ(20, parameters.getGrid_points_per_lambda());

}
TEST_F(Parameters_loadParams_Test, failsForIncorrectValuesInFileOfGoodLength){
    std::string filename = testing_data_dir + "good_length_file_with_params_with_bad_values.txt";
    EXPECT_THROW(parameters.loadParams(filename), std::logic_error);
}
TEST_F(Parameters_loadParams_Test, failsWhenFileDoesNotExist){
    std::string filename = testing_data_dir + "some_non_existing_file.dat";
    EXPECT_THROW(parameters.loadParams(filename), std::ios_base::failure);
}
TEST_F(Parameters_loadParams_Test , failsForFileOfWrongLength){
    std::string filename = testing_data_dir + "4Lines_file_with_params.txt";
    EXPECT_THROW(parameters.loadParams(filename), std::logic_error);
}

//**************loadV()***************
TEST_F(Parameters_loadV_Test, failsIfFileDoesNotExist){
    //try to load some file which does not exist
    EXPECT_THROW(parameters.loadV(testing_data_dir + "some_non_existing_file.dat"), std::ios_base::failure);

}
TEST_F(Parameters_loadV_Test, failsForIncorrectNumberOfRows){
    parameters.nx = 6;
    parameters.nChannels = 1;
    //load file with only 4 values instead of 6
    EXPECT_THROW(parameters.loadV(testing_data_dir + "good_V_one_channel.dat"), std::logic_error);

}
TEST_F(Parameters_loadV_Test, worksGoodForGoodFileOneChannel){
    parameters.nx = 4;
    parameters.nChannels = 1;
    EXPECT_NO_THROW(parameters.loadV(testing_data_dir + "good_V_one_channel.dat"));

    EXPECT_EQ(1.5, std::real(parameters.getVMatrix(0)(0,0)));
    EXPECT_EQ(0, std::imag(parameters.getVMatrix(0)(0,0)));
    EXPECT_EQ(2, std::real(parameters.getVMatrix(1)(0,0)));
    EXPECT_EQ(0, std::imag(parameters.getVMatrix(1)(0,0)));
    EXPECT_EQ(3.5, std::real(parameters.getVMatrix(2)(0,0)));
    EXPECT_EQ(0, std::imag(parameters.getVMatrix(2)(0,0)));
    EXPECT_EQ(4, std::real(parameters.getVMatrix(3)(0,0)));
    EXPECT_EQ(0, std::imag(parameters.getVMatrix(3)(0,0)));



}
TEST_F(Parameters_loadV_Test, worksGoodForGoodFileTwoChannels){
    parameters.nx = 3;
    parameters.nChannels = 2;
    EXPECT_NO_THROW(parameters.loadV(testing_data_dir + "good_V_two_channels_nx=3.txt"));

    EXPECT_DOUBLE_EQ(1, std::real(parameters.getVMatrix(0)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(0)(0,0)));
    EXPECT_DOUBLE_EQ(2, std::real(parameters.getVMatrix(0)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(0)(0,1)));
    EXPECT_DOUBLE_EQ(3, std::real(parameters.getVMatrix(0)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(0)(1,0)));
    EXPECT_DOUBLE_EQ(4, std::real(parameters.getVMatrix(0)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(0)(1,1)));

    EXPECT_DOUBLE_EQ(9.1, std::real(parameters.getVMatrix(2)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(2)(0,0)));
    EXPECT_DOUBLE_EQ(10.1, std::real(parameters.getVMatrix(2)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(2)(0,1)));
    EXPECT_DOUBLE_EQ(11.1, std::real(parameters.getVMatrix(2)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(2)(1,0)));
    EXPECT_DOUBLE_EQ(12.1, std::real(parameters.getVMatrix(2)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getVMatrix(2)(1,1)));


}
//TEST(ParametersInputTest, loadV_failsIfIncorrectNumerOfColumnsInAnyRow);

//***************setXValues()************************
TEST(ParametersInputTest, setXValues_failsIfXMaxLessOrEqualXMin){
    Parameters parameters;
    parameters.xMin = 0;
    parameters.xMax = -1;
    EXPECT_THROW(parameters.setXValues(), std::logic_error);
}
TEST(ParametersInputTest, setXValues_failsIfInvalidDX){
    Parameters parameters;
    parameters.dx = 0;
    EXPECT_THROW(parameters.setXValues(), std::logic_error);
}
TEST(ParametersInputTest, setXValues_failsIfInvalidCombinationOfXMinXMaxDx){
    Parameters parameters;
    parameters.xMin = 0;
    parameters.xMax = 1;
    parameters.dx = 2;
    EXPECT_THROW(parameters.setXValues(), std::logic_error);
}
TEST(ParametersInputTest, setXValues_worksGoodForCorrectValues){
    Parameters parameters;
    parameters.xMin = 0;
    parameters.xMax = 1;
    parameters.dx = 0.5;

    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_EQ(3, parameters.nx);
    EXPECT_DOUBLE_EQ(0, parameters.x(0));
    EXPECT_DOUBLE_EQ(0.5, parameters.x(1));
    EXPECT_DOUBLE_EQ(1, parameters.x(2));
}
//**************x(i)***********************
TEST(ParametersOutputTest, x_worksCorrectForNegativeIndices){
    Parameters parameters;
    parameters.xMin = 0;
    parameters.xMax = 1;
    parameters.dx = 0.5;
    EXPECT_NO_THROW(parameters.setXValues());

    EXPECT_DOUBLE_EQ(1, parameters.x(-1));
    EXPECT_DOUBLE_EQ(0.5, parameters.x(-2));
    EXPECT_DOUBLE_EQ(0, parameters.x(-3));
}

//******************loadE()***********************
class Parameters_loadE_Test: public ::testing::Test{
public:
    Parameters parameters;
};
TEST_F(Parameters_loadE_Test, worksGoodForGoodFileWithDoubleValues){
    parameters.loadE(testing_data_dir + "good_EList_doubles.dat");
    EXPECT_DOUBLE_EQ(1.5, parameters.getE(0));
    EXPECT_DOUBLE_EQ(2.5, parameters.getE(1));
    EXPECT_DOUBLE_EQ(3.41, parameters.getE(2));
}
TEST_F(Parameters_loadE_Test, worksGoodForGoodFileWithIntValues){
    parameters.loadE(testing_data_dir + "good_EList_ints.dat");
    EXPECT_EQ(1, parameters.getE(0));
    EXPECT_EQ(2, parameters.getE(1));
    EXPECT_EQ(3, parameters.getE(2));
};
TEST_F(Parameters_loadE_Test, worksGoodForGoodFileWithMixedDoubleAndIntValues){
    parameters.loadE(testing_data_dir + "good_EList_doubles_ints.dat");
    EXPECT_DOUBLE_EQ(1.5, parameters.getE(0));
    EXPECT_EQ(2, parameters.getE(1));
    EXPECT_DOUBLE_EQ(3.41, parameters.getE(2));
}
TEST_F(Parameters_loadE_Test, failsIfFileDoesNotExist){
    EXPECT_THROW(parameters.loadE(testing_data_dir + "some_non_existing_file.dat"), std::ios_base::failure);
}
TEST_F(Parameters_loadE_Test, failsIfFoundTwoIdenticalValues){
    EXPECT_THROW(parameters.loadE(testing_data_dir + "EList_with_repeats.dat"), std::logic_error);
}
TEST_F(Parameters_loadE_Test, failsIfEmptyFile){
    EXPECT_THROW(parameters.loadE(testing_data_dir + "EList_empty.dat"), std::logic_error);
}

//************loadB()***********************
TEST(ParametersInputTest, loadB_failsIfAnyFileDoesNotExistAndPositiveNSymmetries){
    Parameters parameters;
    parameters.nChannels = 1;
    parameters.nSymmetries = 2;
    EXPECT_THROW(parameters.loadB(testing_data_dir + "some_non_existing_file"), std::logic_error);
}
TEST(ParametersInputTest, loadB_worksGoodForGoodFilesOneChannel){
    Parameters parameters;
    parameters.nChannels = 1;
    parameters.nSymmetries = 2;

    EXPECT_NO_THROW(parameters.loadB(testing_data_dir + "good_one_channel_B"));
    EXPECT_DOUBLE_EQ(1., std::real(parameters.getB(0)(0,0)));
    EXPECT_DOUBLE_EQ(0., std::imag(parameters.getB(0)(0,0)));
    EXPECT_DOUBLE_EQ(-1., std::real(parameters.getB(1)(0,0)));
    EXPECT_DOUBLE_EQ(0., std::imag(parameters.getB(1)(0,0)));

}
TEST(ParametersInputTest, loadB_failsIfAnyFileIsIncorrect){
    Parameters parameters;
    parameters.nChannels = 1;
    parameters.nSymmetries = 3;
    //only B1 and B0 files exist
    EXPECT_THROW(parameters.loadB(testing_data_dir + "good_one_channel_B"), std::logic_error);
}

//************isOpen(nChannel, energy)****************
class Parameters_isOpen_Test: public ::testing::Test{
public:
    Parameters parameters;
    std::string filename;
    void SetUp() override {
        parameters.nx = 3;
        parameters.nChannels = 2;
        filename = testing_data_dir + "good_V_two_channels_nx=3.txt";

    }
};

TEST_F(Parameters_isOpen_Test, trueIfTheChannelIsOpen){
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_TRUE(parameters.isOpen(0, 10));
    EXPECT_TRUE(parameters.isOpen(1, 13));
}
TEST_F(Parameters_isOpen_Test, falseIfTheChannelIsClosed){
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_FALSE(parameters.isOpen(0, 9));
    EXPECT_FALSE(parameters.isOpen(1, 12));
}
TEST_F(Parameters_isOpen_Test, failsIfIncorrectChannelNumber){
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.isOpen(10, 9), std::invalid_argument);
    EXPECT_THROW(parameters.isOpen(-1, 12),std::invalid_argument);
}

//******************kappa(n1, n2, i, E)********************
class Parameters_kappaInt_Test : public ::testing::Test{
public:
    Parameters parameters;
    std::string filename;
    void SetUp() override{
        parameters.xMin = 0;
        parameters.xMax = 1;
        parameters.dx = 0.5;
        parameters.nChannels = 2;
        parameters.unit = 0.5;
        filename = testing_data_dir + "good_V_two_channels_nx=3.txt";
    }
};
TEST_F(Parameters_kappaInt_Test, correctValueFromCorrectArguments){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(2., parameters.kappa(0,0,0,3));
    EXPECT_DOUBLE_EQ(2., parameters.kappa(1,1, 0, 6));
    EXPECT_DOUBLE_EQ(1.2, parameters.kappa(0,0,0,1.72));
}
TEST_F(Parameters_kappaInt_Test, failsIfIncorrectNIndex){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.kappa(5,0,0,3), std::invalid_argument);
    EXPECT_THROW(parameters.kappa(1,-1, 0, 6), std::invalid_argument);
}
TEST_F(Parameters_kappaInt_Test, failsIfIncorrectIIndex){
        EXPECT_NO_THROW(parameters.setXValues());
        EXPECT_NO_THROW(parameters.loadV(filename));
        EXPECT_THROW(parameters.kappa(0,0,-50,3), std::invalid_argument);
        EXPECT_THROW(parameters.kappa(0,0, 10, 6), std::invalid_argument);
}
//***************kappa(n1, n2, x, E)*****************
class Parameters_kappaDouble_Test : public ::testing::Test{
public:
    Parameters parameters;
    std::string filename;
    void SetUp() override{
        parameters.xMin = 0;
        parameters.xMax = 1;
        parameters.dx = 0.5;
        parameters.nChannels = 2;
        parameters.unit = 0.5;
        filename = testing_data_dir + "good_V_two_channels_nx=3.txt";
    }
};
TEST_F(Parameters_kappaDouble_Test, correctValueFromCorrectArguments){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(2., parameters.kappa(0,0,0,3));
}
TEST_F(Parameters_kappaDouble_Test, failsIfIncorrectNIndex){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.kappa(5,0,0.25,3), std::invalid_argument);
    EXPECT_THROW(parameters.kappa(1,-1, 0.25, 6), std::invalid_argument);
}
TEST_F(Parameters_kappaDouble_Test, failsIfIncorrectX){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.kappa(0,0,-50,3), std::invalid_argument);
    EXPECT_THROW(parameters.kappa(0,0, 10, 6), std::invalid_argument);
}

//****************getV(x)******************
class Parameters_getV_Test : public ::testing::Test{
public:
    Parameters parameters;
    std::string filename;
    void SetUp() override{
        parameters.xMin = 0;
        parameters.xMax = 1;
        parameters.dx = 0.5;
        parameters.nChannels = 2;
        filename = testing_data_dir + "good_V_two_channels_nx=3.txt";
    }
};

TEST_F(Parameters_getV_Test, failsForIncorrectX){

    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.getV(-1), std::invalid_argument);
    EXPECT_THROW(parameters.getV(2), std::invalid_argument);
    EXPECT_NO_THROW(parameters.getV(0.3));
}
TEST_F(Parameters_getV_Test, givesCorrectValueForXMin){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));

    EXPECT_DOUBLE_EQ(1, std::real(parameters.getV(0)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0)(0,0)));
    EXPECT_DOUBLE_EQ(2, std::real(parameters.getV(0)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0)(0,1)));
    EXPECT_DOUBLE_EQ(3, std::real(parameters.getV(0)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0)(1,0)));
    EXPECT_DOUBLE_EQ(4, std::real(parameters.getV(0)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0)(1,1)));

}
TEST_F(Parameters_getV_Test, givesCorrectValueForXMax){

    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(9.1, std::real(parameters.getV(1)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(1)(0,0)));
    EXPECT_DOUBLE_EQ(10.1, std::real(parameters.getV(1)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(1)(0,1)));
    EXPECT_DOUBLE_EQ(11.1, std::real(parameters.getV(1)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(1)(1,0)));
    EXPECT_DOUBLE_EQ(12.1, std::real(parameters.getV(1)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(1)(1,1)));
}
TEST_F(Parameters_getV_Test, givesCorrectValueForXOnGrid){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(5, std::real(parameters.getV(0.5)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.5)(0,0)));
    EXPECT_DOUBLE_EQ(6, std::real(parameters.getV(0.5)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.5)(0,1)));
    EXPECT_DOUBLE_EQ(7, std::real(parameters.getV(0.5)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.5)(1,0)));
    EXPECT_DOUBLE_EQ(8, std::real(parameters.getV(0.5)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.5)(1,1)));
}
TEST_F(Parameters_getV_Test, givesCorrectValueForXInterpolated){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(3, std::real(parameters.getV(0.25)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.25)(0,0)));
    EXPECT_DOUBLE_EQ(4, std::real(parameters.getV(0.25)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.25)(0,1)));
    EXPECT_DOUBLE_EQ(5, std::real(parameters.getV(0.25)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.25)(1,0)));
    EXPECT_DOUBLE_EQ(6, std::real(parameters.getV(0.25)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.25)(1,1)));

    EXPECT_DOUBLE_EQ(2, std::real(parameters.getV(0.125)(0,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.125)(0,0)));
    EXPECT_DOUBLE_EQ(3, std::real(parameters.getV(0.125)(0,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.125)(0,1)));
    EXPECT_DOUBLE_EQ(4, std::real(parameters.getV(0.125)(1,0)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.125)(1,0)));
    EXPECT_DOUBLE_EQ(5, std::real(parameters.getV(0.125)(1,1)));
    EXPECT_DOUBLE_EQ(0, std::imag(parameters.getV(0.125)(1,1)));
}

//********************lambda(x, E)*******************************
class Parameters_lambda_Test : public ::testing::Test{
public:
    Parameters parameters;
    std::string filename;
    void SetUp() override{
        parameters.xMin = 0;
        parameters.xMax = 1;
        parameters.dx = 0.5;
        parameters.nChannels = 2;
        parameters.unit = 0.5;
        filename = testing_data_dir + "good_V_two_channels_nx=3.txt";
    }
};
TEST_F(Parameters_lambda_Test, givesCorrectValueForCorrectArguments){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(M_PI/std::sqrt(2.), parameters.lambda(0, 5));
    EXPECT_DOUBLE_EQ(M_PI/std::sqrt(2.), parameters.lambda(0.125, 5.));
    EXPECT_DOUBLE_EQ(2. * M_PI / std::sqrt(14.), parameters.lambda(0.5, 5.));
    EXPECT_DOUBLE_EQ(2. * M_PI / std::sqrt(20.), parameters.lambda(0.125, -5));
}
TEST_F(Parameters_lambda_Test, failsForIncorrectX){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.lambda(-10, 5), std::invalid_argument);
    EXPECT_THROW(parameters.lambda(10, 5), std::invalid_argument);
}

//*************************requiredDx(x, E)*************************
class Parameters_requiredDX_Test : public ::testing::Test{
public:
    Parameters parameters;
    std::string filename;
    void SetUp() override{
        parameters.xMin = 0;
        parameters.xMax = 1;
        parameters.dx = 0.5;
        parameters.nChannels = 2;
        parameters.unit = 1;
        filename = testing_data_dir + "good_V_two_channels_nx=3.txt";
        parameters.grid_points_per_lambda = 1;
    }
};

TEST_F(Parameters_requiredDX_Test, givesCorrectValueForCorrectArgument){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_DOUBLE_EQ(M_PI, parameters.requiredDx(0, 5));

}
TEST_F(Parameters_requiredDX_Test, failsForIncorrectX){
    EXPECT_NO_THROW(parameters.setXValues());
    EXPECT_NO_THROW(parameters.loadV(filename));
    EXPECT_THROW(parameters.requiredDx(-10, 5), std::invalid_argument);
    EXPECT_THROW(parameters.requiredDx(10, 5), std::invalid_argument);
}

//**********************************************************************************************************************
//************************Parameters constructor tests******************************************************************
//TEST(Parameters_constructor_test, worksGoodForCorrectFilesAndParametersWithB){
//    std::vector<std::string> filenamesList = {
//            testing_data_dir + "good_file_with_params.txt",
//            testing_data_dir + "good_EList_doubles.dat",
//            testing_data_dir + "good_V_two_channels_nx=3.txt",
//            testing_data_dir + "good_B_two_channels"
//    };
//    Parameters p(filenamesList);
//}

//TEST(Parameters_constructor_test, worksGoodForCorrectFilesAndParametersWithoutB);
//TEST(Parameters_constructor_test, exitsIfIncorrectSizeOfListOfFilenames){
//        std::vector<std::string> filenamesList = {
//            testing_data_dir + "good_file_with_params.txt",
//            testing_data_dir + "good_EList_doubles.dat",
//    };
//    EXPECT_EXIT(Parameters p(filenamesList), ::testing::ExitedWithCode(-2), "");
//}
//TEST(Parameters_constructor_test, exitsIfIncorrectValuesInAnyOfDataFiles);
//TEST(Parameters_constructor_test, exitsIfCannotReadAnyOfTheDataFiles);


#endif //NUMEROV_PARAMETERSTESTS_H
