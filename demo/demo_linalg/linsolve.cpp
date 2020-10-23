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
 | ( === ) |   SYNOPSIS   : Solve a linear system A*X=B                   |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+==============+" << std::endl;
    std::cout << "|   LINSOLVE   |" << std::endl;
    std::cout << "+==============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;

    // FLOAT
    auto Af = rand<float>(m);
    auto Bf = rand<float>(m,2);
    matrix<float> Xf = linsolve(Af,Bf);
    disp( norm(mtimes(Af,Xf)-Bf,"inf") );
    
    Af = cat(1,Af,Af); // vertical concatenation
    Bf = cat(1,Bf,Bf);
    Xf = linsolve(Af,Bf);
    disp( norm(mtimes(Af,Xf)-Bf,"inf") );
    
    Af = rand<float>(m,n);
    Bf = rand<float>(m,2);
    Xf = linsolve(Af,Bf);
    disp( norm(mtimes(Af,Xf)-Bf,"inf") );
    
    // DOUBLE
    auto Ad = rand(m);
    auto Bd = rand(m,2);
    matrix<> Xd = linsolve(Ad,Bd);
    disp( norm(mtimes(Ad,Xd)-Bd,"inf") );
    
    Ad = cat(1,Ad,Ad);
    Bd = cat(1,Bd,Bd);
    Xd = linsolve(Ad,Bd);
    disp( norm(mtimes(Ad,Xd)-Bd,"inf") );
    
    Ad = rand(m,n);
    Bd = rand(m,2);
    Xd = linsolve(Ad,Bd);
    disp( norm(mtimes(Ad,Xd)-Bd,"inf") );
    
    // COMPLEX FLOAT
    auto Ac = rand<float>(m) + M_1If*rand<float>(m);
    auto Bc = rand<float>(m,2) + M_1If*rand<float>(m,2);
    auto Xc = linsolve(Ac,Bc);
    disp( norm(mtimes(Ac,Xc)-Bc,"inf") );
    
    Ac = cat(1,Ac,Ac);
    Bc = cat(1,Bc,Bc);
    Xc = linsolve(Ac,Bc);
    disp( norm(mtimes(Ac,Xc)-Bc,"inf") );
    
    Ac = rand<float>(m,n) + M_1If*rand<float>(m,n);
    Bc = rand<float>(m,2) + M_1If*rand<float>(m,2);
    Xc = linsolve(Ac,Bc);
    disp( norm(mtimes(Ac,Xc)-Bc,"inf") );
    
    // COMPLEX DOUBLE
    auto Az = rand(m) + M_1I*rand(m);
    auto Bz = rand(m,2) + M_1I*rand(m,2);
    auto Xz = linsolve(Az,Bz);
    disp( norm(mtimes(Az,Xz)-Bz,"inf") );
    
    Az = cat(1,Az,Az);
    Bz = cat(1,Bz,Bz);
    Xz = linsolve(Az,Bz);
    disp( norm(mtimes(Az,Xz)-Bz,"inf") );

    Az = rand(m,n) + M_1I*rand(m,n);
    Bz = rand(m,2) + M_1I*rand(m,2);
    Xz = linsolve(Az,Bz);
    disp( norm(mtimes(Az,Xz)-Bz,"inf") );
    
    disp("done !");
    return 0;
}
