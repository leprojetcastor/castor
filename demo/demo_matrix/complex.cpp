/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : complex.cpp                                   |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Use and functions of complex values           |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+===============+");
    disp("|    COMPLEX    |");
    disp("+===============+");
    
    matrix<> Re = {1,0,-1,0};
    matrix<> Im = {0,1,0,-1};

    disp("Imaginary number matrix ('1i') : ");
    disp(M_1I);

    disp("Complex matrix : ");
    auto Ac = Re + M_1I*Im;
    disp(Ac);
    
    disp("Amplitude of a complex matrix : ");
    disp(abs(Ac));
    
    disp("Phase of a complex matrix (in radian) : ");
    disp(angle(Ac));
    
    disp("Conjugate matrix : ");
    disp(conj(Ac));
    
    disp("Exponential of a complex matrix : ");
    disp(exp(M_1I*Re*M_PI));
    
    disp("Real part and immaginary part : ");
    disp(real(Ac));
    disp(imag(Ac));
    
    disp("Square roots :");
    disp(sqrt(Ac));
    
    disp("Pow :");
    disp(pow(Ac,2));
    
    disp("Sum :");
    disp(sum(Ac));
    
    disp("Prod :");
    disp(prod(Ac));
    
    disp("Norm :");
    disp(norm(Ac));
    
    
    disp("done !");
    return 0;
}
