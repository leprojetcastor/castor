/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : qrsvd.cpp                                     |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute the rank-revealing SVD                |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|     QRSVD     |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // FLOAT
    auto Af = rand(m,k), Bf = rand(k,n);
    matrix<float> Uf, Vf;
    std::tie(Uf,Vf) = qrsvd(Af,Bf,1e-3);
    disp(norm(mtimes(Af,Bf)-mtimes(Uf,Vf),"inf"));

    // DOUBLE
    auto Ad = rand(m,k), Bd = rand(k,n);
    matrix<double> Ud, Vd;
    std::tie(Ud,Vd) = qrsvd(Ad,Bd,1e-3);
    disp(norm(mtimes(Ad,Bd)-mtimes(Ud,Vd),"inf"));

    // COMPLEX FLOAT
    auto Ac = rand<float>(m,k) + M_1If*rand<float>(m,k);
    auto Bc = rand<float>(k,n) + M_1If*rand<float>(k,n);
    matrix<std::complex<float>> Uc, Vc;
    std::tie(Uc,Vc) = qrsvd(Ac,Bc,1e-3);
    disp(norm(mtimes(Ac,Bc)-mtimes(Uc,Vc),"inf"));

    // COMPLEX DOUBLE
    auto Az = rand(m,k) + M_1I*rand(m,k);
    auto Bz = rand(k,n) + M_1I*rand(k,n);
    matrix<std::complex<double>> Uz, Vz;
    std::tie(Uz,Vz) = qrsvd(Az,Bz,1e-3);
    disp(norm(mtimes(Az,Bz)-mtimes(Uz,Vz),"inf"));

    disp("done !");
    return 0;
}
