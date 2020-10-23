/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : manipulation.cpp                              |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Manipulates matrices and related data         |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+===================+");
    disp("|   MANIPULATIONS   |");
    disp("+===================+");

    matrix<> V =  {0,1,2,M_PI};
    matrix<> A = {{0,1,2,M_PI},
        {4,5,6,7},
        {8,9,10,11}};
    auto B = A;
    
    disp("Cast matrix values from type double to int :");
    disp(cast<int>(A));
    disp(cast(A,int()));

    disp("Vertical concatenation of scalars to build matrix :");
    disp(cat(1,0,0));
    disp(vertcat(0,0));
    disp("Horizontal concatenation of scalars to build matrix :");
    disp(cat(2,0,0));
    disp(horzcat(0,0));
    
    disp("Vertical concatenation of matrix with scalar :");
    disp(cat(1,0,ones(3,1)));
    disp(vertcat(0,ones(3,1)));
    disp("Horizontal concatenation of matrix with scalar :");
    disp(cat(2,0,ones(1,3)));
    disp(horzcat(0,ones(1,3)));
    
    disp("Vertical concatenation of matrices :");
    disp(cat(1,A,eye(size(A))));
    disp(vertcat(A,eye(size(A))));
    disp("Horizontal concatenation of matrices :");
    disp(cat(2,A,eye(size(A))));
    disp(horzcat(A,eye(size(A))));
    
    disp("Clear matrix to free memory :");
    clear(B);
    disp(B);

    disp("Find linear indices of non-zero elements :");
    disp(find(eye(3,4)));
    
    disp("All linear indices of the matrix :");
    disp(all(A));
    
    disp("Column indices of the matrix :");
    disp(col(A));

    disp("Row indices of the matrix :");
    disp(row(A));
    
    disp("Get matrix values from linear indices :");
    disp(get(A,all(A)));

    disp("Set matrix values from linear indices :");
    B = A;
    set(B,all(B),M_PI);
    disp(B);
    set(B,all(B),rand(1,12));
    disp(B);

    disp("Get matrix values from bi-linear indices :");
    disp(get(A,row(A),0));
    disp(get(A,0,col(A)));
    disp(get(A,range(0,2),range(0,2)));
    
    disp("Set matrix values from bi-linear indices :");
    set(B,row(B),0,-1);
    disp(B);
    set(B,0,col(B),-1);
    disp(B);
    set(B,row(B),col(B),rand(3,4));
    disp(B);
    
    disp("Convert linear to bi-linear indices :");
    matrix<std::size_t> idx,jdx;
    std::tie(idx,jdx) = ind2sub({2,3},range(0,6));
    disp(idx);
    disp(jdx);
    
    disp("Convert bi-linear to linear indices :");
    disp(sub2ind({2,3},{0,1,2},{0,1,2}));

    disp("Test if matrix is empty :");
    disp(isempty(matrix<>()));
    disp(isempty(A));
    
    disp("Test if matrix are equal :");
    disp(isequal(A,A));
    disp(isequal(A,A+1));

    disp("Test if matrix is mono-dimensional :");
    disp(isvector(ones(1,3)));
    disp(isvector(ones(3,1)));
    disp(isvector(ones(3,3)));
    
    disp("Reshape matrix form conserving row-major layout :");
    disp(reshape(A,4,3));

    disp("Transpose matrix :");
    disp(transpose(A));
    
    disp("Resize matrix :");
    disp(resize(A,4,5));
    
    disp("Values of the matrix :");
    disp(values(A));
    
    disp("done !");
    return 0;
}
