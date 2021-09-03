/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : builder.cpp                                   |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Matrix initializers and builders              |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    disp("+=====================+");
    disp("|      BUILDERS       |");
    disp("+=====================+");

    disp("Build scalar logical matrix :");
    matrix<logical> s = (logical)true;
    disp(s);
    
    disp("Build row int vector :");
    matrix<int> X = {0,1,2,3};
    disp(X);
    
    disp("Build float matrix :");
    matrix<float> A = {{0,1,2,M_PI},
        {4,5,6,7},
        {8,9,10,11}};
    disp(A);
    
    disp("Copy float to default (double) matrix :");
    matrix<> B = A;
    disp(B);
    
    disp("Identity matrix of logical with same size as A:");
    disp(eye<logical>(size(A)));
    disp("Identity matrix of double:");
    disp(eye(4,3));
    
    disp("Ones matrix of logical with same size as A:");
    disp(ones<logical>(size(A)));
    disp("Ones matrix of double:");
    disp(ones(4,3));

    disp("Uniform random matrix of float with same size as A:");
    disp(rand<float>(size(A)));
    disp("Uniform random matrix of double:");
    disp(rand(4,3));
    
    disp("Normal random matrix of float with same size as A:");
    disp(randn<float>(size(A)));
    disp("Normal random matrix of double:");
    disp(randn(4,3));

    disp("Random permutation:");
    disp(randperm(5));
    disp(randperm(5,3));
    disp(randperm(5,3,true));
    
    disp("Range :");
    disp(range(0,4));

    disp("Zeros matrix of int with same size as A:");
    disp(rand<int>(size(A)));
    disp("Zeros matrix of double:");
    disp(zeros(4,3));
    
    disp("Generate row vector with 1-step values from 0 to 5 :");
    disp(colon(0,5));
    disp("Generate row vector with 1.5-step values from 0 to 7 :");
    disp(colon(0,1.5,7));
    
    disp("Linear repartion of 5 values, from 1 to 0 :");
    disp(linspace(1,0,5));
    disp("Logarithmic repartion of 5 values, from 10^0 to 10^3 :");
    disp(logspace(0,3,4));
    
    disp("Generate diagonal matrix from vector :");
    disp(diag(X));
    disp("Upper diagonal :");
    disp(diag(X,1));
    disp("Lower diagonal :");
    disp(diag(X,-2));

    disp("Extract diagonal from matrix :");
    disp(diag(A));
    disp("Upper diagonal :");
    disp(diag(A,1));
    disp("Lower diagonal :");
    disp(diag(A,-2));
    
    disp("Hann window");
    disp(hann(5));

    disp("done !");
    return 0;
}
