/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : eig.cpp                                       |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Find eigenvalues and eigenvectors             |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|      EIG      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // FLOAT
    auto Af = 1-eye<float>(m);
    matrix<float> Ef, Vf, Uf;
    std::tie(Ef,Uf) = eig(Af,"left");
    Uf = transpose(Uf);
    disp(norm(mtimes(Uf,Af)-mtimes(diag(Ef),Uf),"inf"));
    std::tie(Ef,Vf) = eig(Af,"right");
    disp(norm(mtimes(Af,Vf)-mtimes(Vf,diag(Ef)),"inf"));

    // DOUBLE
    auto Ad = 1-eye(m);
    matrix<> Ed, Vd, Ud;
    std::tie(Ed,Ud) = eig(Ad,"left");
    Ud = transpose(Ud);
    disp(norm(mtimes(Ud,Ad)-mtimes(diag(Ed),Ud),"inf"));
    std::tie(Ed,Vd) = eig(Ad,"right");
    disp(norm(mtimes(Ad,Vd)-mtimes(Vd,diag(Ed)),"inf"));

    // COMPLEX FLOAT
    auto Ac = rand<float>(m) + M_1If*rand<float>(m);
    matrix<std::complex<float>> Ec, Vc, Uc;
    std::tie(Ec,Uc) = eig(Ac,"left");
    Uc = conj(transpose(Uc));
    disp(norm(mtimes(Uc,Ac)-mtimes(diag(Ec),Uc),"inf"));
    std::tie(Ec,Vc) = eig(Ac,"right");
    disp(norm(mtimes(Ac,Vc)-mtimes(Vc,diag(Ec)),"inf"));

    // COMPLEX DOUBLE
    auto Az = rand(m) + M_1I*rand(m);
    matrix<std::complex<double>> Ez, Vz, Uz;
    std::tie(Ez,Uz) = eig(Az,"left");
    Uz = conj(transpose(Uz));
    disp(norm(mtimes(Uz,Az)-mtimes(diag(Ez),Uz),"inf"));
    std::tie(Ez,Vz) = eig(Az,"right");
    disp(norm(mtimes(Az,Vz)-mtimes(Vz,diag(Ez)),"inf"));

    disp("done !");
    return 0;
}
