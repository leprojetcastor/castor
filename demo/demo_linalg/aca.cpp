/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : aca.cpp                                       |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Adaptive cross approximation                  |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|      ACA      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;

    // DOUBLE
    matrix<> Ad, Bd;
    matrix<> Md = mtimes(rand(m,k),rand(k,n));
    std::tie(Ad,Bd) = aca(Md,1e-3);
    disp(size(Ad));
    disp(norm(Md-mtimes(Ad,Bd),"inf"));

    // COMPLEX DOUBLE
    matrix<std::complex<double>> Mz = mtimes(rand(m,k),rand(k,n))*M_1I;
    matrix<std::complex<double>> Az, Bz;
    std::tie(Az,Bz) = aca(Mz,1e-3);
    disp(size(Az));
    disp(norm(Mz-mtimes(Az,Bz),"inf"));

    disp("done !");
    return 0;
}
