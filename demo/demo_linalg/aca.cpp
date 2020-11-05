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
    // Parameters
    std::size_t m=100, n=200, k=10;
    bool flag;
    
    std::cout << "+============+" << std::endl;
    std::cout << "|   DOUBLE   |" << std::endl;
    std::cout << "+============+" << std::endl;

    // BLOC 1x1
    matrix<> Ad, Bd, Md;
    Md = mtimes(rand(m,k),rand(k,n));
    std::tie(Ad,Bd,flag) = aca(Md);
    disp(size(Ad));
    disp(rank(Md));
    disp(norm(Md-mtimes(Ad,Bd),"inf"));
    
    // BLOC 2x2
    matrix<> Zd = zeros(m,n);
    Md = vertcat(horzcat(Md,Zd),horzcat(Zd,Md));
    std::tie(Ad,Bd,flag) = aca(Md);
    disp(size(Ad));
    disp(rank(Md));
    disp(norm(Md-mtimes(Ad,Bd),"inf"));
    
    std::cout << "+============+" << std::endl;
    std::cout << "|   COMPLEX  |" << std::endl;
    std::cout << "+============+" << std::endl;

    // BLOC 1x1
    matrix<std::complex<double>> Az, Bz, Mz;
    Mz = mtimes(rand(m,k),rand(k,n)) + M_1I*mtimes(rand(m,k),rand(k,n));
    std::tie(Az,Bz,flag) = aca(Mz);
    disp(size(Az));
    disp(rank(Mz));
    disp(norm(Mz-mtimes(Az,Bz),"inf"));
    
    // BLOC 2x2
    matrix<std::complex<double>> Zz = zeros(m,n);
    Mz = vertcat(horzcat(Mz,Zz),horzcat(Zz,Mz));
    std::tie(Az,Bz,flag) = aca(Mz);
    disp(size(Az));
    disp(rank(Mz));
    disp(norm(Mz-mtimes(Az,Bz),"inf"));
    
    std::cout << "+=================+" << std::endl;
    std::cout << "|  RECOMPRESION   |" << std::endl;
    std::cout << "+=================+" << std::endl;
    
    // BLOC 1x1
    Md = diag(logspace(0,-14,m));
    std::tie(Ad,Bd,flag) = aca(Md);
    disp(size(Ad));
    disp(flag);
    disp(norm(Md-mtimes(Ad,Bd),"inf"));
    
    // Fix tolerance
    std::tie(Ad,Bd,flag) = aca(Ad,Bd,1e-3);
    disp(size(Ad));
    disp(flag);
    disp(norm(Md-mtimes(Ad,Bd),"inf"));
    
    // Fix rank
    std::tie(Ad,Bd,flag) = aca(Ad,Bd,1e-3,10);
    disp(size(Ad));
    disp(flag);    
    disp(norm(Md-mtimes(Ad,Bd),"inf"));

    disp("done !");
    return 0;
}
