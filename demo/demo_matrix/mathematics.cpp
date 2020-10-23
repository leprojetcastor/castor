/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : mathematics.cpp                               |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Basic mathematical function of matrix         |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    disp("+=====================+");
    disp("|     MATHEMATICS     |");
    disp("+=====================+");
    
    matrix<> V = {-M_PI,-2,-1,0,1,2,M_PI};
    std::cout << "original : " << V << std::endl;
    std::cout << "abs      : " << abs(V) << std::endl;
    std::cout << "acos     : " << acos(V) << std::endl;
    std::cout << "acosd    : " << acosd(V) << std::endl;
    std::cout << "acosh    : " << acosh(V) << std::endl;
    std::cout << "asin     : " << asin(V) << std::endl;
    std::cout << "asind    : " << asind(V) << std::endl;
    std::cout << "asinh    : " << asinh(V) << std::endl;
    std::cout << "atan     : " << atan(V) << std::endl;
    std::cout << "atand    : " << atand(V) << std::endl;
    std::cout << "atanh    : " << atanh(V) << std::endl;
    std::cout << "ceil     : " << ceil(V) << std::endl;
    std::cout << "cos      : " << cos(V) << std::endl;
    std::cout << "cosd     : " << cosd(V*180/M_PI) << std::endl;
    std::cout << "cosh     : " << cosh(V) << std::endl;
    std::cout << "deg2rad  : " << deg2rad(V*180/M_PI) << std::endl;
    std::cout << "exp      : " << exp(V) << std::endl;
    std::cout << "floor    : " << floor(V) << std::endl;
    std::cout << "log      : " << log(V) << std::endl;
    std::cout << "log2     : " << log2(V) << std::endl;
    std::cout << "log10    : " << log10(V) << std::endl;
    std::cout << "pow      : " << pow(V,V) << std::endl;
    std::cout << "pow      : " << pow(V,2) << std::endl;
    std::cout << "pow      : " << pow(2,V) << std::endl;
    std::cout << "rad2deg  : " << rad2deg(V) << std::endl;
    std::cout << "round    : " << round(V) << std::endl;
    std::cout << "sign     : " << sign(V) << std::endl;
    std::cout << "sin      : " << sin(V) << std::endl;
    std::cout << "sind     : " << sind(V*180/M_PI) << std::endl;
    std::cout << "sinh     : " << sinh(V) << std::endl;
    std::cout << "sqrt     : " << sqrt(V) << std::endl;
    std::cout << "tan      : " << tan(V) << std::endl;
    std::cout << "tand     : " << tand(V*180/M_PI) << std::endl;
    std::cout << "tanh     : " << tanh(V) << std::endl;
    
    disp("done !");
    return 0;
}
