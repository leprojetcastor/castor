/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : svd.cpp                                       |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute the Singular Value Decomposition      |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|      SVD      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200;

    // FLOAT
    auto Af = rand<float>(m,n);
    matrix<float> Sf, Uf, Vf;
    std::tie(Sf,Uf,Vf) = svd(Af,"vect");
    disp(norm(Af-mtimes(Uf,mtimes(diag(Sf),Vf)),"inf"));

    Af = rand<float>(n,m);
    std::tie(Sf,Uf,Vf) = svd(Af,"vect");
    disp(norm(Af-mtimes(Uf,mtimes(diag(Sf),Vf)),"inf"));

    // DOUBLE
    auto Ad = rand(m,n);
    matrix<> Sd, Ud, Vd;
    std::tie(Sd,Ud,Vd) = svd(Ad,"vect");
    disp(norm(Ad-mtimes(Ud,mtimes(diag(Sd),Vd)),"inf"));

    Ad = rand(n,m);
    std::tie(Sd,Ud,Vd) = svd(Ad,"vect");
    disp(norm(Ad-mtimes(Ud,mtimes(diag(Sd),Vd)),"inf"));

    // COMPLEX FLOAT
    auto Ac = rand<float>(m,n)*M_1If;
    matrix<std::complex<float>> Uc, Vc;
    std::tie(Sf,Uc,Vc) = svd(Ac,"vect");
    disp(norm(Ac-mtimes(Uc,mtimes(diag(Sf),Vc)),"inf"));

    Ac = rand<float>(n,m)*M_1If;
    std::tie(Sf,Uc,Vc) = svd(Ac,"vect");
    disp(norm(Ac-mtimes(Uc,mtimes(diag(Sf),Vc)),"inf"));

    // COMPLEX DOUBLE
    auto Az = rand(m,n)*M_1I;
    matrix<std::complex<double>> Uz, Vz;
    std::tie(Sd,Uz,Vz) = svd(Az,"vect");
    disp(norm(Az-mtimes(Uz,mtimes(diag(Sd),Vz)),"inf"));

    Az = rand(n,m)*M_1I;
    std::tie(Sd,Uz,Vz) = svd(Az,"vect");
    disp(norm(Az-mtimes(Uz,mtimes(diag(Sd),Vz)),"inf"));

    disp("done !");
    return 0;
}
