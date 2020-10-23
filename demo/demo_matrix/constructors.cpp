/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : constructors.cpp                              |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Examples of matrix class constructors         |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    disp("+=====================+");
    disp("|    CONSTRUCTORS     |");
    disp("+=====================+");

    disp("Empty logical matrix : ");
    matrix<logical> A1;
    disp(A1);

    disp("Scalar int matrix from fundamental type constant or variable :");
    matrix<int> A2(3);
    disp(A2);
    
    disp("Row float matrix from initializer list :");
    matrix<float> A3({0,1,2,M_PI});
    disp(A3);

    disp("Double Matrix from bi-initializer list :");
    matrix<double> A4({{0,1,2,M_PI},{4,5,6,7},{8,9,10,11}});
    disp(A4);

    disp("Complex Matrix from value filling:");
    matrix<std::complex<double>> A5(2,3,std::complex<double>(1,M_PI));
    disp(A5);

    disp("Row double Matrix from std::vector of int :");
    matrix<double> A6(std::vector<int>({0, 1, 2, 3}));
    disp(A6);
    
    disp("Default (double) Matrix from size and initializer list, with auto cast :");
    matrix<> A7(2,2,{0, 1, 2, 3});
    disp(A7);
    matrix<> A8(2,2,{{0, 1},{2,3}});
    disp(A8);
    
    disp("Default (double) Matrix from matrix copy with auto cast :");
    matrix<> A9(A8);
    disp(A9);
    
    disp("done !");
    return 0;
}
