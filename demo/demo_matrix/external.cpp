/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : external.cpp                                  |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : External operators of matrix class            |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+======================+");
    disp("|  EXTERNAL OPERATORS  |");
    disp("+======================+");
    
    matrix<int>     I = {1,1,0,0};
    matrix<logical> J = {0,1,0,1};
    matrix<float>   A = eye(3,4);
    matrix<>        V = {0,1,2,3};
    
    disp("operator+");
    std::cout << "A      : " << V << std::endl;
    std::cout << "A+1    : " << (V + 1) << std::endl;
    std::cout << "1+A    : " << (1 + V) << std::endl;
    std::cout << "A+v    : " << (V + sum(V)) << std::endl;
    std::cout << "v+A    : " << (sum(V) + V) << std::endl;
    std::cout << "A+B<>  : " << (V + matrix<int>(1,4,1)) << std::endl;
    std::cout << "A<>+B  : " << (matrix<int>(1,4,1) + V) << std::endl;
    std::cout << "[]+1   : " << (matrix<>() + 1) << std::endl;
    std::cout << "1+[]   : " << (1 + matrix<>()) << std::endl;
    
    disp("operator-");
    std::cout << "A      : " << V << std::endl;
    std::cout << "-A     : " << (-V) << std::endl;
    std::cout << "A-1    : " << (V - 1) << std::endl;
    std::cout << "1-A    : " << (1 - V) << std::endl;
    std::cout << "A-B<>  : " << (V - matrix<int>(1,4,1)) << std::endl;
    std::cout << "A<>-B  : " << (matrix<int>(1,4,1) - V) << std::endl;
    
    disp("operator*");
    std::cout << "A      : " << V << std::endl;
    std::cout << "A*2    : " << (V * 2) << std::endl;
    std::cout << "2*A    : " << (2 * V) << std::endl;
    std::cout << "A*B<>  : " << (V * matrix<int>(1,4,2)) << std::endl;
    std::cout << "A<>*B  : " << (matrix<int>(1,4,2) * V) << std::endl;

    disp("operator/");
    std::cout << "A      : " << V << std::endl;
    std::cout << "A/2    : " << (V / 2) << std::endl;
    std::cout << "2/A    : " << (2 / V) << std::endl;
    std::cout << "A/B<>  : " << (V / matrix<int>(1,4,2)) << std::endl;
    std::cout << "A<>/B  : " << (matrix<int>(1,4,2) / V) << std::endl;
    
    disp("operator!");
    std::cout << "A      : " << I << std::endl;
    std::cout << "!A     : " << (!I) << std::endl;

    disp("operator&&");
    std::cout << "A      : " << I << std::endl;
    std::cout << "B      : " << J << std::endl;
    std::cout << "A&&B   : " << (I&&J) << std::endl;
    std::cout << "A&&v   : " << (I&&0) << std::endl;
    std::cout << "v&&B   : " << (0&&J) << std::endl;

    disp("operator||");
    std::cout << "A      : " << I << std::endl;
    std::cout << "B      : " << J << std::endl;
    std::cout << "A||B   : " << (I||J) << std::endl;
    std::cout << "A||1   : " << (I||1) << std::endl;
    std::cout << "1||B   : " << (1||J) << std::endl;
    
    disp("operator@=");
    std::cout << "A      : " << I << std::endl;
    std::cout << "B      : " << J << std::endl;
    std::cout << "A==B   : " << (I==J) << std::endl;
    std::cout << "A==1   : " << (I==1) << std::endl;
    std::cout << "1==A   : " << (1==I) << std::endl;
    std::cout << "A!=B   : " << (I!=J) << std::endl;
    std::cout << "A!=1   : " << (I!=1) << std::endl;
    std::cout << "1!=A   : " << (1!=I) << std::endl;
    std::cout << "A<=B   : " << (I<=J) << std::endl;
    std::cout << "A<=0   : " << (I<=0) << std::endl;
    std::cout << "1<=A   : " << (1<=I) << std::endl;
    std::cout << "A<B    : " << (I<J) << std::endl;
    std::cout << "A<1    : " << (I<1) << std::endl;
    std::cout << "0<A    : " << (0<I) << std::endl;
    std::cout << "A>=B   : " << (I>=J) << std::endl;
    std::cout << "A>=1   : " << (I>=1) << std::endl;
    std::cout << "0>=A   : " << (0>=I) << std::endl;
    std::cout << "A>B    : " << (I>J) << std::endl;
    std::cout << "A>0    : " << (I>0) << std::endl;
    std::cout << "1>A    : " << (1>I) << std::endl;
    
    disp("Dimensions of a matrix :");
    std::cout << "length(A) : " << length(A) << std::endl;
    std::cout << "nnz(A)    : " << nnz(A) << std::endl;
    std::cout << "numel(A)  : " << numel(A) << std::endl;
    std::cout << "size(A,0) : " << size(A,0) << std::endl;
    std::cout << "size(A,1) : " << size(A,1) << std::endl;
    std::cout << "size(A,2) : " << size(A,2) << std::endl;
    std::cout << "size(A)   : " << size(A) << std::endl;
    
    disp("done !");
    return 0;
}
