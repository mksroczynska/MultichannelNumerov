//
// Created by martas on 06.12.17.
//

#include <string>
#include "Parameters.h"

/*!
 *\file
 * \brief Definitions of Parameters class methods.
 */

Parameters::Parameters(std::vector<std::string> filenames) {
    if (filenames.size() == 3 or filenames.size() == 4) {
        try {
            loadParams(filenames.at(0));
        } catch (std::ios_base::failure &exception) {
            std::cout << "Error loading values of parameters:\n" << exception.what();
            exit(-2);
        } catch (std::logic_error &exception) {
            std::cout << "Error loading values of parameters:\n" << exception.what();
            exit(-1);
        }
        try{
            setXValues();
        }catch(std::logic_error &exception){
            std::cout << "Error creating grid:\n" << exception.what();
        }
        try {
            loadE(filenames.at(1));
        } catch (std::ios_base::failure &exception) {
            std::cout << "Error loading values of energies:\n" << exception.what();
            exit(-2);
        } catch (std::logic_error &exception) {
            std::cout << "Error loading values of energies:\n" << exception.what();
            exit(-1);
        }
        try {
            loadV(filenames.at(2));
        } catch (std::ios_base::failure &exception) {
            std::cout << "Error loading values of potential:\n" << exception.what();
            exit(-2);
        } catch (std::logic_error &exception) {
            std::cout << "Error loading values of potential:\n" << exception.what();
            exit(-1);
        }
    } else {
        std::cout << "Wrong size of the list of filenames: " << filenames.size() << ", required: 3 or 4";
        exit(-2);
    }
    if (nSymmetries > 0) {
        if (filenames.size() == 4) {
            try {
                loadB(filenames.at(3));
            } catch (std::logic_error &exception) {
                std::cout << "Error loading values of B:\n" << exception.what();
                exit(-1);
            }
        } else {
            std::cout << "No required filename for files with B\n";
            exit(-2);
        }

    }

}

arma::cx_mat Parameters::getVMatrix(int i) const {

    if(i < -nx or i > nx)
        throw std::invalid_argument("getVMatrix: wrong index");
    if(i >= 0)
        return V.slice(i);
    else
        return V.slice(nx + i);
}

double Parameters::getE(int i) const {
    return EList.at(i);
}

int Parameters::NX() const {
    return nx;
}

void Parameters::loadParams(std::string filename) {
    std::ifstream fileP(filename);
    if (fileP.good()) {
        if (checkNumberOfRowsInFile(filename, 7) or checkNumberOfRowsInFile(filename, 6) or
            checkNumberOfRowsInFile(filename, 5)) {
            fileP >> xMin;
            fileP >> xMax;
            if (xMax <= xMin)
                throw std::logic_error("Invalid values of xMin and xMax: "
                                       + std::to_string(xMin) + ", " + std::to_string(xMax));
            fileP >> dx;
            if (dx <= 0)
                throw std::logic_error("Invalid value of dx: " + std::to_string(dx));
            if (xMax - xMin < dx)
                throw std::logic_error("Invalid values: xMax - xMin > dx");
            fileP >> unit;
            if (unit <= 0)
                throw std::logic_error("Invalid value of unit: " + std::to_string(unit));
            fileP >> nChannels;
            if (nChannels <= 0)
                throw std::logic_error("Invalid number of channels: " + std::to_string(nChannels));
            nSymmetries = 0;
            if (checkNumberOfRowsInFile(filename, 6) or checkNumberOfRowsInFile(filename, 7)) {
                fileP >> nSymmetries;
                if (nSymmetries < 0)
                    throw std::logic_error("Invalid number of symmetries: " + std::to_string(nSymmetries));
            }
            if (checkNumberOfRowsInFile(filename, 7)) {
                fileP >> grid_points_per_lambda;
                if (grid_points_per_lambda <= 0 and nSymmetries > 0)
                    throw std::logic_error("Invalid number of grid points per deBroglie wavelength: " +
                                           std::to_string(grid_points_per_lambda));
            }
        } else {
            throw std::logic_error("File with parameters has wrong length");
        }
        fileP.close();
    } else {
        throw std::ios_base::failure("No required file with parameters: " + filename);
    }

    std::cout << "***************************************************\n";
    std::cout << "Loaded values of the parameters:\n";
    std::cout << "xMin: " << xMin << "\nxMax: " << xMax << "\ndx: " << dx <<
              "\nunit: " << unit << "\nnumber of channels: " << nChannels <<
              "\nnumber of symmetries: " << nSymmetries << '\n' << "grid points per lambda: " << grid_points_per_lambda;
    std::cout << "\n***************************************************\n";
}

void Parameters::loadE(std::string filename) {

    std::ifstream fileE(filename);
    if (fileE.good()) {
        double ev;
        std::string line;
        while (std::getline(fileE, line)) {
            ev = std::stod(line);
            auto it = std::find(std::begin(EList), std::end(EList), ev);
            if (it != std::end(EList)) {
                throw std::logic_error(
                        "Incorrect list of energies: more than one occurence of value " + std::to_string(ev));
            }
            EList.push_back(ev);
        }
        nE = EList.size();
        if (nE < 1) throw std::logic_error("Empty file with energies");
        fileE.close();
    } else {
        throw std::ios_base::failure("No required file with energies: " + filename);
    }
}

void Parameters::loadV(std::string filename) {

    std::ifstream fileV(filename);

    if (fileV.good()) {
        if (checkNumberOfRowsInFile(filename, NX())) {
            arma::cx_cube V_c(nChannels, nChannels, NX());
            arma::cx_mat V_j(nChannels, nChannels);

            for (auto j = 0; j < NX(); ++j) {
                for (int i = 0; i < nChannels; ++i) {
                    for (int k = 0; k < nChannels; ++k) {
                        double vik;
                        fileV >> vik;
                        V_j(i, k) = vik;
                    }
                }
                V_c.slice(j) = V_j;
            }
            V = V_c;
        } else {
            throw std::logic_error("File with V has wrong number of rows ");
        }
    } else {
        throw std::ios_base::failure("No required file with potential: " + filename);
    }
}


bool Parameters::isOpen(int nChannel, double energy) const {
    if (nChannel >= nChannels or nChannel < 0)
        throw std::invalid_argument("Cannot check if channel " + std::to_string(nChannel) + "is open: " +
                                    "wrong channel index");
    return (energy - ((getVMatrix(nx - 1))(nChannel, nChannel)).real() > 0);
}

double Parameters::kappa(int n1, int n2, int i, double E) const {
    if (n1 >= nChannels or n1 < 0 or n2 >= nChannels or n2 < 0)
        throw std::invalid_argument(
                "Cannot calculate kappa for (" + std::to_string(n1) + "," + std::to_string(n2) + "): wrong index");
    if (i >= nx)
        throw std::invalid_argument("Cannot calculate kappa for x(i): wrong value of i");
    auto v = ((getVMatrix(i))(n1, n2)).real();
    if (E - v > 0)
        return std::sqrt(1./unit * (E - v));
    else
        return std::sqrt(1./unit * (v - E));
}

void Parameters::loadB(std::string filename) {
    arma::cx_cube B_temp(nChannels, nChannels, nSymmetries);

    for (int i = 0; i < nSymmetries; ++i) {
        auto filenameB = filename + std::to_string(i) + ".dat";
        bool success_loading_B = (B_temp.slice(i)).load(filenameB);
        if (not success_loading_B)
            throw std::logic_error("File with B_" + std::to_string(i) + " matrix does not exist");
    }

    B = B_temp;

}

arma::cx_mat Parameters::getV(double x) const {
    // 1) x out of range
    if (x < xMin or x > xMax)
        throw std::invalid_argument("Interpolation: incorrect value of argument: " + std::to_string(x));

    // 2) x on grid
    auto it = std::find(xValues.begin(), xValues.end(), x);
    if (it != xValues.end()) {
        return getVMatrix(it - xValues.begin());
    }
        // 3) interpolate
    else {
        auto index_2 =
                std::find_if(xValues.begin(), xValues.end(), [x](double value) { return x < value; }) - xValues.begin();
        auto index_1 = index_2 - 1;
        arma::cx_mat Y1 = V.slice(index_1);
        arma::cx_mat Y2 = V.slice(index_2);
        auto x1 = xValues.at(index_1);
        auto x2 = xValues.at(index_2);
        arma::cx_mat A = (Y2 - Y1) / (x2 - x1);
        arma::cx_mat B = Y2 - x2 * A;

        return A * x + B;
    }
}

double Parameters::lambda(double x, double E) const {
    if (x < xMin or x > xMax)
        throw std::invalid_argument("Cannot calculate lambda for x =" + std::to_string(x) + ": wrong value of x");
    arma::cx_mat Id(nChannels, nChannels, arma::fill::eye);
    arma::mat absEV = arma::abs(E * Id - getV(x));
    double max_element_EV = absEV.max();
    return 2. * M_PI / std::sqrt(1./unit * max_element_EV);
}

double Parameters::kappa(int n1, int n2, double x, double E) const {
    if (n1 >= nChannels or n1 < 0 or n2 >= nChannels or n2 < 0)
        throw std::invalid_argument(
                "Cannot calculate kappa for (" + std::to_string(n1) + "," + std::to_string(n2) + "): wrong index");
    if (x < xMin or x > xMax)
        throw std::invalid_argument("Cannot calculate kappa for x =" + std::to_string(x) + ": wrong value of x");
    auto v = ((getV(x))(n1, n2)).real();
    if (E - v > 0)
        return std::sqrt(1./unit  * (E - v));
    else
        return std::sqrt(1./unit * (v - E));
}

double Parameters::requiredDx(double x, double E) const {
    return lambda(x, E) / grid_points_per_lambda;
}


void Parameters::setXValues() {
    if (xMax <= xMin)
        throw std::logic_error("Invalid values of xMin and xMax");
    if (dx <= 0)
        throw std::logic_error("Invalid value of dx");
    if (xMax - xMin < dx)
        throw std::logic_error("Invalid values of xMin, xMax and dx");
    for (auto value = xMin; value <= xMax; value += dx) {
        xValues.push_back(value);
    }
    nx = xValues.size();
}

bool Parameters::checkNumberOfRowsInFile(std::string filename, const int required_number_of_rows) {
    std::ifstream file(filename);
    int calculated_number_of_rows = 0;
    std::string line;
    while (std::getline(file, line))
        calculated_number_of_rows++;
    return calculated_number_of_rows == required_number_of_rows;
}































