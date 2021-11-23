/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : getting_started.cpp                           |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Running first demo for matrix library         |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    disp("+=====================+");
    disp("|   GETTING STARTED   |");
    disp("+=====================+");

    // Configure documnentation file
    documentationFiles =
    {
        "/usr/local/include/castor/matrix.hpp"  // Full name of your file matrix.hpp
    };

    // Documentation of matrix class
    help("matrix");

    // Timer
    tic();
    while (toc(0)<0.314159) {}
    toc();

    // Define matrix by values
    matrix<> A = {{ 1.0,  2.0,  3.0,  4.0},
                  { 5.0,  6.0,  7.0,  8.0},
                  { 9.0, 10.0, 11.0, 12.0}};

    // Define matrix using builder
    auto B = eye(3,4);

    // Standard operation
    auto C = A + B;

    // Display
    disp("C = A+B :");
    disp(C);

    disp("done !");
    return 0;
}
