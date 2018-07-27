//
// Created by martas on 06.12.17.
//

#ifndef NOMICROMOTION_PARAMETERS_H
#define NOMICROMOTION_PARAMETERS_H


#include <cmath>
#include <vector>
#include <string>
#include <exception>
#include <gtest/gtest_prod.h>
#include "armadillo"

/*!
* \file
 * \brief Definition of Parameters class
 *
 * This file contains the definition of Parameters class.
*/

class Parameters {

    /*!
     * \brief The smallest value of x.
     */
    double xMin;
    /*!
     * \brief The largest value of x.
     */
    double xMax;
    /*!
     * \brief Step: x_{i+1} = x_{i} + dx.
     */
    double dx;
    /*!
     * \brief The list of x values from xMin to xMax
     */
    std::vector<double> xValues;
    /*!
     * \brief List of the energies for which the calculations can be performed within this set of parameters.
     */
    std::vector<double> EList;
    /*!
     * \brief hbar^2/2m_a.
     */
    double unit;
    /*!
     * \brief Number of channels.
     */
    int nChannels;
    /*!
     * \brief Number of symmetries.
     */
    int nSymmetries;
    /*!
     * \brief Cube storing V matrices for each x
     */
    arma::cx_cube V;
    /*!
     * \brief EList length.
     */
    unsigned int nE;
    /*
     * \brief How many grid points.
     */
    int nx;
    /*!
     * \brief Cube storing B matrix for each symmetry.
     */
    arma::cx_cube B;
    /*!
     * \brief How many grid points of Numerov iteration should be per de Broglie wavelength .
     */
    int grid_points_per_lambda;
    /*!
     * \brief Reading the values of parameters from the file generated in Mathematica.
     */
public:
    void loadParams(std::string = "Params.txt");

    /*!
     * \brief Setting xValues
     */
    void setXValues();
    FRIEND_TEST(ParametersInputTest, setXValues_failsIfXMaxLessOrEqualXMin);
    FRIEND_TEST(ParametersInputTest, setXValues_failsIfInvalidDX);
    FRIEND_TEST(ParametersInputTest, setXValues_failsIfInvalidCombinationOfXMinXMaxDx);
    FRIEND_TEST(ParametersInputTest, setXValues_worksGoodForCorrectValues);
    /*!
     * \brief Reading the values of energies from the file generated in Mathematica.
     */
    void loadE(std::string filename = "E.dat");

    /*!
     * \brief Reading the values of V from the file generated in Mathematica.
     */
    void loadV(std::string filename = "V.dat");
    friend class Parameters_loadV_Test;
    FRIEND_TEST(Parameters_loadV_Test, failsForIncorrectNumberOfRows);
    FRIEND_TEST(Parameters_loadV_Test, worksGoodForGoodFileOneChannel);
    FRIEND_TEST(Parameters_loadV_Test, worksGoodForGoodFileTwoChannels);

    /*!
     * \brief Reading the values of B from the file generated in Mathematica.
     */
    void loadB(std::string filename = "B");
    FRIEND_TEST(ParametersInputTest, loadB_failsIfAnyFileDoesNotExistAndPositiveNSymmetries);
    FRIEND_TEST(ParametersInputTest, loadB_worksGoodForGoodFilesOneChannel);
    FRIEND_TEST(ParametersInputTest, loadB_failsIfAnyFileIsIncorrect);

    bool checkNumberOfRowsInFile(std::string filename, const int required_number_of_columns);

public:
    Parameters() = default;
    ~Parameters() = default;
    /*!
     * \brief From a given directory takes all the needed values and creates Parameters object.
     */
    explicit Parameters(std::vector<std::string> filenames);

    /*!
     * \brief V matrix for a given x_i.
     */
    arma::cx_mat getVMatrix(int) const;

    double getE(int) const;


    int NX() const;


    double getXMin() const {
        return xMin;
    }

    double getXMax() const {
        return xMax;
    }

    double getDx() const {
        return dx;
    }

    double x(int i) const {
        if(i >= 0)
            return xValues.at(i);
        else
            return xValues.at(nx + i);
    }
    FRIEND_TEST(ParametersOutputTest, x_worksCorrectForNegativeIndices);


    double getUnit() const {
        return unit;
    }

    int getNChannels() const {
        return nChannels;
    }

    int getNE() const {
        return nE;
    }

    arma::cx_mat getB(int i) const {
        if(i >= nSymmetries and i > 0)
            throw std::invalid_argument("B_" + std::to_string(i) + " does not exist since the number of symmetries = " + std::to_string(nSymmetries));
        if(i == 0 and nSymmetries == 0) return arma::cx_mat(nChannels, nChannels, arma::fill::zeros);
        return B.slice(i);
    }


    int getNSymmetries() const {
        return nSymmetries;
    }

    int getGrid_points_per_lambda() const
    {
        return grid_points_per_lambda;
    }

    arma::cx_mat Id() const{
        return arma::cx_mat(nChannels, nChannels, arma::fill::eye);
    }

    /*!
     * \brief Check if the channel is open.
     */
    bool isOpen(int nChannel, double energy) const;
    friend class Parameters_isOpen_Test;

    double kappa(int n1, int n2, int i, double E) const;
    friend class Parameters_kappaInt_Test;

    double kappa(int n1, int n2, double x, double E) const;
    friend class Parameters_kappaDouble_Test;

    /*!
     * \brief Linear interpolation of V (works also for V given on non-constant grid if needed)
     */
    arma::cx_mat getV(double x) const;
    friend class Parameters_getV_Test;

    /*!
     * \brief de Broglie length for a given potential and x
     */
    double lambda(double x, double E) const;
    friend class Parameters_lambda_Test;


    double requiredDx(double x, double E) const;
    friend class Parameters_requiredDX_Test;

};


#endif //NOMICROMOTION_PARAMETERS_H
