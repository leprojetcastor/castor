/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : rank.cpp                                      |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Compute rank of matrix                        |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/linalg.hpp"
#include "castor/matrix.hpp"

using namespace castor;

int main()
{
    std::cout << "+===============+" << std::endl;
    std::cout << "|      RANK     |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // Matrix size
    std::size_t m=100, n=200, k=20;
    
    // FLOAT
    auto Af = mtimes(rand<float>(m,k),rand<float>(k,n));
    disp(rank(Af));

    // DOUBLE
    auto Ad = mtimes(rand(m,k),rand(k,n));
    disp(rank(Ad));

    // COMPLEX FLOAT
    auto Ac = mtimes(rand<float>(m,k),rand<float>(k,n))*M_1If;
    disp(rank(Ac));

    // COMPLEX DOUBLE
    auto Az = mtimes(rand(m,k),rand(k,n))*M_1I;
    disp(rank(Az));

    disp("done !");
    return 0;
}
