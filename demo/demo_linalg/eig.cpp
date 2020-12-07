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
    matrix<std::complex<float>> Ec, Uc, Vc;
    std::tie(Ec,Uc) = eig(Af,"left");
    Uc = conj(transpose(Uc));
    disp(norm(mtimes(Uc,Af)-mtimes(diag(Ec),Uc),"inf"));
    std::tie(Ec,Vc) = eig(Af,"right");
    disp(norm(mtimes(Af,Vc)-mtimes(Vc,diag(Ec)),"inf"));

    // DOUBLE
    auto Ad = 1-eye(m);
    matrix<std::complex<double>> Ez, Uz, Vz;
    std::tie(Ez,Uz) = eig(Ad,"left");
    Uz = conj(transpose(Uz));
    disp(norm(mtimes(Uz,Ad)-mtimes(diag(Ez),Uz),"inf"));
    std::tie(Ez,Vz) = eig(Ad,"right");
    disp(norm(mtimes(Ad,Vz)-mtimes(Vz,diag(Ez)),"inf"));
    
    // COMPLEX FLOAT
    auto Ac = rand<float>(m) + M_1If*rand<float>(m);
    std::tie(Ec,Uc) = eig(Ac,"left");
    Uc = conj(transpose(Uc));
    disp(norm(mtimes(Uc,Ac)-mtimes(diag(Ec),Uc),"inf"));
    std::tie(Ec,Vc) = eig(Ac,"right");
    disp(norm(mtimes(Ac,Vc)-mtimes(Vc,diag(Ec)),"inf"));

    // COMPLEX DOUBLE
    auto Az = rand(m) + M_1I*rand(m);
    std::tie(Ez,Uz) = eig(Az,"left");
    Uz = conj(transpose(Uz));
    disp(norm(mtimes(Uz,Az)-mtimes(diag(Ez),Uz),"inf"));
    std::tie(Ez,Vz) = eig(Az,"right");
    disp(norm(mtimes(Az,Vz)-mtimes(Vz,diag(Ez)),"inf"));

    // FLOAT - GENERALIZED
    Af      = rand<float>(m);
    auto Bf = rand<float>(m);
    std::tie(Ec,Uc) = eig(Af,Bf,"left");
    Uc = conj(transpose(Uc));
    disp(norm(mtimes(Uc,Af) - mtimes(diag(Ec),mtimes(Uc,Bf)),"inf"));
    std::tie(Ec,Vc) = eig(Af,Bf,"right");
    disp(norm(mtimes(Af,Vc) - mtimes(Bf,mtimes(Vc,diag(Ec))),"inf"));

    // // DOUBLE - GENERALIZED
    Ad      = rand<double>(m);
    auto Bd = rand<double>(m);
    std::tie(Ez,Uz) = eig(Ad,Bd,"left");
    Uz = conj(transpose(Uz));
    disp(norm(mtimes(Uz,Ad) - mtimes(diag(Ez),mtimes(Uz,Bd)),"inf"));
    std::tie(Ez,Vz) = eig(Ad,Bd,"right");
    disp(norm(mtimes(Ad,Vz) - mtimes(Bd,mtimes(Vz,diag(Ez))),"inf"));

    // COMPLEX FLOAT - GENERALIZED
    Ac      = rand<float>(m) + M_1If*rand<float>(m);
    auto Bc = rand<float>(m) + M_1If*rand<float>(m);
    std::tie(Ec,Uc) = eig(Ac,Bc,"left");
    Uc = conj(transpose(Uc));
    disp(norm(mtimes(Uc,Ac) - mtimes(diag(Ec),mtimes(Uc,Bc)),"inf"));
    std::tie(Ec,Vc) = eig(Ac,Bc,"right");
    disp(norm(mtimes(Ac,Vc) - mtimes(Bc,mtimes(Vc,diag(Ec))),"inf"));
    
    // COMPLEX DOUBLE - GENERALIZED
    Az      = rand<double>(m) + M_1I*rand<double>(m);
    auto Bz = rand<double>(m) + M_1I*rand<double>(m);
    std::tie(Ez,Uz) = eig(Az,Bz,"left");
    Uz = conj(transpose(Uz));
    disp(norm(mtimes(Uz,Az) - mtimes(diag(Ez),mtimes(Uz,Bz)),"inf"));
    std::tie(Ez,Vz) = eig(Az,Bz,"right");
    disp(norm(mtimes(Az,Vz) - mtimes(Bz,mtimes(Vz,diag(Ez))),"inf"));

    disp("done !");
    return 0;
}
