#include <iostream>
#include <vector>
#include <fstream>
#include "Parameters.h"
#include "ConstantGridSolver.h"
#include "NonconstantGridSolver.h"

using namespace std;

std::vector<std::string> read_file(std::string filename) {
    std::ifstream file(filename);
    std::vector<std::string> result;
    std::string line;
    while (std::getline(file, line)) {
        result.push_back(line);
    }
    return result;
}

int main(int argc, char **argv) {

    std::string directory = argv[1];//"//home//martas//DoktoratoweData//Chain//NumerovAE//R=0.5//";
            std::cout << argv[1];
    //std::string directory = "//home//martas//DoktoratoweData//AIcollisions//noMicromotion//AI//1//";
    //std::string directory = "//home//martas//Dropbox//Doktoratowe//NumerovAE//R=0.2//";
    auto params_data = directory + "Params.dat";
    auto E_data = directory + "E.dat";
    auto v_filenames = read_file(directory + "names.txt");
    int i = 1;
    for (auto v_data: v_filenames) {
        auto filenames = {
                params_data,
                E_data,
                directory + v_data,
                directory + "B"
        };
        Parameters parameters(filenames);
        ConstantGridSolver solver(parameters);


        solver.solveForEnergies(directory + "S//S_" + std::to_string(i) + ".dat");
        // solver.solveForEnergies(directory+ "S//S.dat");
        i++;
    }


}
