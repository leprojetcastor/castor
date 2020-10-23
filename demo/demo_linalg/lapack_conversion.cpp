/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : linsolve.cpp                                  |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Conversion between the matrix format and the  |
 |  `---'  |                Fortran column-ordering                       |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===================+" << std::endl;
    std::cout << "| CONVERSION LAPACK |" << std::endl;
    std::cout << "+===================+" << std::endl;

    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // Conversion between matrices of various floating point precision
    matrix<double> A     = rand(m,n);
    std::vector<float> V = mat2lpk<float>(A,numel(A));
    matrix<double> B     = lpk2mat<double>(V,size(A,1),size(A,2));
    std::cout << "|A - B|_inf = " << norm(A-B,"inf") << std::endl << std::endl;

    // Conversion from the std::complex  to the
    // C/Fortran compatible complex type

    // testing single precision
    clpk cl = std::complex<float>(1,1); // cl = (1,1)
    cl.r    = 2; // cl = (2,1)
    cl.i    = 3; // cl = (2,3)
    std::cout << "cl = " << std::flush;
    disp(cl);
    std::complex<float> cc = cl; // conversion to the C++ type
    std::cout << "cc = ";
    disp(cc);
    std::cout  << std::endl;

    // same but for double precision
    zlpk zl = std::complex<double>(1,1);
    zl.r = 4;
    zl.i = 5;
    std::cout << "zl = " << std::flush;
    disp(zl);
    std::complex<double> zc = zl;
    std::cout << "zl = " << std::flush;
    disp(zc);

    disp("done !");
    return 0;
}
