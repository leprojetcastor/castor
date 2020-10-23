/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : view.cpp                                      |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : View matrices and assign values               |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+====================+");
    disp("|        VIEW        |");
    disp("+====================+");
    
    matrix<> A = {{0,1,2,M_PI},
        {4,5,6,7},
        {8,9,10,11}};
    matrix<> B;
    matrix<> V = {-1,-1,-1,-1};
    
    disp("Original matrix :");
    disp(A);
    
    
    disp("+-----------------+");
    disp("| LINEAR INDEXING |");
    disp("+-----------------+");

    disp("Evaluate A(L) with initialization list :");
    B = eval(A({1,3,5}));
    disp(B);
    
    disp("Evaluate A(L) with range :");
    B = eval(A(range(0,4)));
    disp(B);
    
    disp("Evaluate A(L) with all indices (linear form) :");
    B = eval(A(all(A)));
    disp(B);

    disp("Assign A(L) with a scalar value :");
    B = A;
    B(range(0,4)) = 0;
    disp(B);
    
    disp("Assign A(L) with an initialization list :");
    B = A;
    B(range(0,4)) = {1,1,1,1};
    disp(B);

    disp("Assign A(L) with a matrix :");
    B = A;
    B(range(0,4)) = eye(1,4);
    disp(B);
    
    
    disp("+--------------------+");
    disp("| BI-LINEAR INDEXING |");
    disp("+--------------------+");

    disp("Evaluate A(I,J) with initialization lists :");
    B = eval(A({1,2},{1,2,3}));
    disp(B);

    disp("Evaluate A(I,J) with range :");
    B = eval(A(range(0,3),range(0,3)));
    disp(B);

    disp("Evaluate A(I,J) with all indices (bi-linear form):");
    B = eval(A(row(A),col(A)));
    disp(B);
    
    disp("Assign A(I,J) with a scalar value :");
    B = A;
    B({0,1},col(B)) = 0;
    disp(B);
    
    disp("Assign A(I,J) with an initialization list :");
    B = A;
    B({0,1},col(B)) = {{1,1,1,1},{1,1,1,1}};
    disp(B);

    disp("Assign A(I,J) with a matrix :");
    B = A;
    B({0,1},col(B)) = eye(2,4);
    disp(B);

    
    disp("done !");
    return 0;
}
