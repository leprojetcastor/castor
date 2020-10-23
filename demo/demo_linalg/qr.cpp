/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : qr.cpp                                        |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute the QR decomposition                  |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      QR       |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // Matrix size
    std::size_t m=100, n=200, k=20;

    // FLOAT
    auto Af = rand<float>(m,n);
    matrix<float> Qf, Rf;
    std::tie(Qf,Rf) = qr(Af);
    disp(norm(Af-mtimes(Qf,Rf),"inf"));

    Af = rand<float>(n,m);
    std::tie(Qf,Rf) = qr(Af);
    disp(norm(Af-mtimes(Qf,Rf),"inf"));

    // DOUBLE
    auto Ad = rand(m,n);
    matrix<double> Qd, Rd;
    std::tie(Qd,Rd) = qr(Ad);
    disp(norm(Ad-mtimes(Qd,Rd),"inf"));

    Ad = rand(n,m);
    std::tie(Qd,Rd) = qr(Ad);
    disp(norm(Ad-mtimes(Qd,Rd),"inf"));

    // COMPLEX FLOAT
    auto Ac = rand<float>(m,n)*M_1If;
    matrix<std::complex<float>> Qc, Rc;
    std::tie(Qc,Rc) = qr(Ac);
    disp(norm(Ac-mtimes(Qc,Rc),"inf"));

    Ac = rand<float>(n,m)*M_1If;
    std::tie(Qc,Rc) = qr(Ac);
    disp(norm(Ac-mtimes(Qc,Rc),"inf"));

    // COMPLEX DOUBLE
    auto Az = rand(m,n)*M_1I;
    matrix<std::complex<double>> Qz, Rz;
    std::tie(Qz,Rz) = qr(Az);
    disp(norm(Az-mtimes(Qz,Rz),"inf"));

    Az = rand(n,m)*M_1I;
    std::tie(Qz,Rz) = qr(Az);
    disp(norm(Az-mtimes(Qz,Rz),"inf"));

    disp("done !");
    return 0;
}
