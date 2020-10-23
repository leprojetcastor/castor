/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : lu.cpp                                        |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute the Lower/Upper factorization         |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|       LU      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // FLOAT
    matrix<float> Lf, Uf;
    auto Af = rand<float>(m,n);
    std::tie(Lf,Uf) = lu(Af);
    disp( norm(mtimes(Lf,Uf)-Af,"inf") );
    
    Af = rand<float>(n,m);
    std::tie(Lf,Uf) = lu(Af);
    disp( norm(mtimes(Lf,Uf)-Af,"inf") );
    
    // DOUBLE
    matrix<> Ld, Ud;
    auto Ad = rand(m,n);
    std::tie(Ld,Ud) = lu(Ad);
    disp( norm(mtimes(Ld,Ud)-Ad,"inf") );
    
    Ad = rand(n,m);
    std::tie(Ld,Ud) = lu(Ad);
    disp( norm(mtimes(Ld,Ud)-Ad,"inf") );
    
    // COMPLEX FLOAT
    matrix<std::complex<float>> Lc, Uc;
    auto Ac = rand<float>(m,n) + M_1If*rand<float>(m,n);
    std::tie(Lc,Uc) = lu(Ac);
    disp( norm(mtimes(Lc,Uc)-Ac,"inf") );
    
    Ac = rand<float>(n,m) + M_1If*rand<float>(n,m);
    std::tie(Lc,Uc) = lu(Ac);
    disp( norm(mtimes(Lc,Uc)-Ac,"inf") );
    
    // COMPLEX DOUBLE
    matrix<std::complex<double>> Lz, Uz;
    auto Az = rand(m,n) + M_1I*rand(m,n);
    std::tie(Lz,Uz) = lu(Az);
    disp( norm(mtimes(Lz,Uz)-Az,"inf") );

    Az = rand(n,m) + M_1I*rand(n,m);
    std::tie(Lz,Uz) = lu(Az);
    disp( norm(mtimes(Lz,Uz)-Az,"inf") );

    disp("done !");
    return 0;
}
