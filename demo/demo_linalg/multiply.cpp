/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : multiply.cpp                                  |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute the matrix-matrix product             |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|   MULTIPLY    |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // FLOAT
    matrix<float> Af = rand<float>(1,m); // pick random matrices
    matrix<float> Bf = rand<float>(m,1);
    float Cf = 0;
    for (std::size_t l=0; l<m; ++l) {Cf+=Af(l)*Bf(l);}
    std::cout << "(float)   |A*B - C| = " << std::flush;
    disp( norm(mtimes(Af,Bf)-Cf,"inf") ); // multiply and compare to ref. C
    std::cout << std::endl;
    
    // DOUBLE
    matrix<> Ad = rand(1,m);
    matrix<> Bd = rand(m,1);
    double Cd = 0;
    for (std::size_t l=0; l<m; ++l) {Cd+=Ad(l)*Bd(l);}
    std::cout << "(double)  |A*B - C| = " << std::flush;
    disp( norm(mtimes(Ad,Bd)-Cd,"inf") );
    std::cout << std::endl;
    
    // COMPLEX FLOAT
    auto Ac = rand<float>(1,m) + M_1If*rand<float>(1,m);
    auto Bc = rand<float>(m,1) + M_1If*rand<float>(m,1);
    std::complex<float> Cc = 0;
    for (std::size_t l=0; l<m; ++l) {Cc+=Ac(l)*Bc(l);}
    std::cout << "(cfloat)  |A*B - C| = " << std::flush;
    disp( norm(mtimes(Ac,Bc)-Cc,"inf") );
    std::cout << std::endl;

    // COMPLEX DOUBLE
    auto Az = rand(1,m) + M_1I*rand(1,m);
    auto Bz = rand(m,1) + M_1I*rand(m,1);
    std::complex<double> Cz = 0;
    for (std::size_t l=0; l<m; ++l) {Cz+=Az(l)*Bz(l);}
    std::cout << "(cdouble) |A*B - C| = " << std::flush;
    disp( norm(mtimes(Az,Bz)-Cz,"inf") );
    std::cout << std::endl;

    disp("done !");
    return 0;
}
