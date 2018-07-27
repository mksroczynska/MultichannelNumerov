#include <iostream>
#include <vector>
#include <fstream>
#include "Parameters.h"
#include "ConstantGridSolver.h"

//#include "NonconstantGridSolver.h"
using namespace std;

int main() {
    std::string testing_data_dir = "//home//martas//Dropbox//Doktoratowe//AIcollisions//MultichannelNumerovAIWithTests_1//MultichannelNumerovAIWithTests//TestData//";

    std::vector<std::string> filenamesList = {
            testing_data_dir + "good_file_with_params.txt",
            testing_data_dir + "good_EList_doubles.dat",
            testing_data_dir + "good_V_two_channels_nx=3.txt",
            testing_data_dir + "good_B_two_channels"
    };
    Parameters p(filenamesList);
}
