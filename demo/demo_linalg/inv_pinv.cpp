/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : inv_pinv.cpp                                  |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute inverse or pseudo-inverse of matrix   |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|      INV      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // FLOAT
    auto Af = rand<float>(m);
    auto Xf = inv(Af);
    disp( norm(mtimes(Af,Xf)-eye<float>(m),"inf") );
    
    // DOUBLE
    auto Ad = rand(m);
    auto Xd = inv(Ad);
    disp( norm(mtimes(Ad,Xd)-eye(m),"inf") );
    
    // COMPLEX FLOAT
    auto Ac = rand<float>(m) + M_1If*rand<float>(m);
    auto Xc = inv(Ac);
    disp( norm(mtimes(Ac,Xc)-eye<float>(m),"inf") );
    
    // COMPLEX DOUBLE
    auto Az = rand(m) + M_1I*rand(m);
    auto Xz = inv(Az);
    disp( norm(mtimes(Az,Xz)-eye(m),"inf") );
    
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|     PINV      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // FLOAT
    Af = rand<float>(m,n);
    Xf = pinv(Af);
    disp(norm(mtimes(mtimes(Xf,Af),Xf)-Xf,"inf"));
    disp(norm(mtimes(mtimes(Af,Xf),Af)-Af,"inf"));
    
    Af = rand<float>(n,m);
    Xf = pinv(Af);
    disp(norm(mtimes(mtimes(Xf,Af),Xf)-Xf,"inf"));
    disp(norm(mtimes(mtimes(Af,Xf),Af)-Af,"inf"));
    
    // DOUBLE
    Ad = rand(m,n);
    Xd = pinv(Ad);
    disp(norm(mtimes(mtimes(Xd,Ad),Xd)-Xd,"inf"));
    disp(norm(mtimes(mtimes(Ad,Xd),Ad)-Ad,"inf"));
    
    Ad = rand(n,m);
    Xd = pinv(Ad);
    disp(norm(mtimes(mtimes(Xd,Ad),Xd)-Xd,"inf"));
    disp(norm(mtimes(mtimes(Ad,Xd),Ad)-Ad,"inf"));
    
    // COMPLEX FLOAT
    Ac = rand<float>(m,n) + M_1If*rand<float>(m,n);
    Xc = pinv(Ac);
    disp(norm(mtimes(mtimes(Xc,Ac),Xc)-Xc,"inf"));
    disp(norm(mtimes(mtimes(Ac,Xc),Ac)-Ac,"inf"));
    
    Ac = rand<float>(n,m) + M_1If*rand<float>(n,m);
    Xc = pinv(Ac);
    disp(norm(mtimes(mtimes(Xc,Ac),Xc)-Xc,"inf"));
    disp(norm(mtimes(mtimes(Ac,Xc),Ac)-Ac,"inf"));
    
    // COMPLEX DOUBLE
    Az = rand(m,n) + M_1I*rand(m,n);
    Xz = pinv(Az);
    disp(norm(mtimes(mtimes(Xz,Az),Xz)-Xz,"inf"));
    disp(norm(mtimes(mtimes(Az,Xz),Az)-Az,"inf"));
    
    Az = rand(n,m) + M_1I*rand(n,m);
    Xz = pinv(Az);
    disp(norm(mtimes(mtimes(Xz,Az),Xz)-Xz,"inf"));
    disp(norm(mtimes(mtimes(Az,Xz),Az)-Az,"inf"));
    
    disp("done !");
    return 0;
}
