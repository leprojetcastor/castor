/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : algorithms.cpp                                |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Standard algorithms for matrix class          |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+==================+");
    disp("|    ALGORITHMS    |");
    disp("+==================+");
    
    matrix<> A = {{0,1,2,M_PI},
        {4,5,6,7},
        {8,9,10,11}};
    A = A + 1;
    
    disp("Original matrix :");
    disp(A);
    
    
    disp("+---------------+");
    disp("|    MINIMUM    |");
    disp("+---------------+");
    
    disp("Minimum value :");
    disp(min(A));

    disp("Index of minimum value (linear indexing) :");
    disp(argmin(A));
    
    disp("Minimum values for each column :");
    disp(min(A,1));
    
    disp("Indices of minimum values for each column (bi-linear indexing) :");
    disp(argmin(A,1));
    
    disp("Minimum values for each row : ");
    disp(min(A,2));
    
    disp("Indices of minimum values for each row (bi-linear indexing) :");
    disp(argmin(A,2));

    disp("Minimum values between two independant matrix :");
    disp(minimum(A,A+1));
    
    disp("Minimum values between matrix and scalar :");
    disp(minimum(A,5.5));
    disp("--------------------");
    disp(minimum(5,A));

    
    disp("+---------------+");
    disp("|    MAXIMUM    |");
    disp("+---------------+");
    
    disp("Maximum value :");
    disp(max(A));

    disp("Index of maximum value (linear indexing) :");
    disp(argmax(A));
    
    disp("Maximum values for each column :");
    disp(max(A,1));
    
    disp("Indices of maximum values for each column (bi-linear indexing) :");
    disp(argmax(A,1));
    
    disp("Maximum values for each row :");
    disp(max(A,2));
    
    disp("Indices of maximum values for each row (bi-linear indexing) :");
    disp(argmax(A,2));

    disp("Maximum values between two independant matrix :");
    disp(maximum(A,A+1));
    
    disp("Maximum values between matrix and scalar :");
    disp(maximum(A,5.5));
    disp("--------------------");
    disp(maximum(5,A));
    
    
    disp("+---------------+");
    disp("|      SORT     |");
    disp("+---------------+");
    
    disp("Sorted data in ascending order :");
    disp(sort(eye(3,4)));
    
    disp("Indices of sorted data in ascending order (linear indexing) :");
    disp(argsort(eye(3,4)));

    disp("Sorted data in ascending order for each column :");
    disp(sort(eye(3,4),1));

    disp("Indices of sorted data in ascending order for each column (bi-linear indexing) :");
    disp(argsort(eye(3,4),1));
    
    disp("Sorted data in ascending order for each row :");
    disp(sort(eye(3,4),2));
    
    disp("Indices of sorted data in ascending order for each row (bi-linear indexing) :");
    disp(argsort(eye(3,4),2));
    

    disp("+------------------+");
    disp("|        SUM       |");
    disp("+------------------+");

    disp("Sum of each element :");
    disp(sum(eye(3,4)));
    
    disp("Sum of each element for each column :");
    disp(sum(eye(3,4),1));
    
    disp("Sum of each element for each row :");
    disp(sum(eye(3,4),2));
    
    
    disp("+------------------+");
    disp("|      CUMSUM      |");
    disp("+------------------+");

    disp("Cumulative sum of each element :");
    disp(cumsum(eye(3,4)));
    
    disp("Cumulative sum of each element for each column :");
    disp(cumsum(eye(3,4),1));
    
    disp("Cumulative sum of each element for each row :");
    disp(cumsum(eye(3,4),2));
    
    
    disp("+------------------+");
    disp("|      PRODUCT     |");
    disp("+------------------+");

    disp("Product of each element :");
    disp(prod(2*ones(3,4)));
    
    disp("Product of each element for each column :");
    disp(prod(2*ones(3,4),1));
    
    disp("Product of each element for each row :");
    disp(prod(2*ones(3,4),2));
    
    
    disp("+------------------+");
    disp("|    CUMPRODUCT    |");
    disp("+------------------+");

    disp("Cumulative product of each element :");
    disp(cumprod(2*ones(3,4)));
    
    disp("Cumulative product of each element for each column :");
    disp(cumprod(2*ones(3,4),1));
    
    disp("Cumulative product of each element for each row :");
    disp(cumprod(2*ones(3,4),2));

    
    disp("+------------------+");
    disp("|    DIFFERENCE    |");
    disp("+------------------+");

    disp("Difference of each element :");
    disp(diff(A));
    
    disp("Difference of each element for each column :");
    disp(diff(A,1));
    
    disp("Difference of each element for each row :");
    disp(diff(A,2));

    
    disp("+------------------+");
    disp("|       MEAN       |");
    disp("+------------------+");

    disp("Mean value for all element :");
    disp(mean(eye(3,4)));
    
    disp("Mean value of each column :");
    disp(mean(eye(3,4),1));
    
    disp("Mean value of each row :");
    disp(mean(eye(3,4),2));

    
    disp("+------------------+");
    disp("|      MEDIAN      |");
    disp("+------------------+");
    
    disp("Median value for all element :");
    disp(median(A));
    
    disp("Median value of each column :");
    disp(median(A,1));
    
    disp("Median value of each row :");
    disp(median(A,2));

    
    disp("+--------------------+");
    disp("|      VARIANCE      |");
    disp("+--------------------+");
    
    disp("Variance of all element :");
    disp(variance(A));
    
    disp("Variance of each column :");
    disp(variance(A,1));
    
    disp("Variance of each row :");
    disp(variance(A,2));
    
    
    disp("+--------------------+");
    disp("| STANDARD DEVIATION |");
    disp("+--------------------+");

    disp("Standard deviation for all element :");
    disp(stddev(A));

    disp("Standard deviation of each column :");
    disp(stddev(A,1));

    disp("Standard deviation of each row :");
    disp(stddev(A,2));
    
    
    disp("+--------------------+");
    disp("|        NORMS       |");
    disp("+--------------------+");

    disp("1-Norm of matrix all values :");
    disp(norm(A,"1"));

    disp("2-Norm of matrix all values :");
    disp(norm(A,"2"));

    disp("Infinite norm of all matrix values :");
    disp(norm(A,"inf"));
    
    disp("Infinite norm of matrix values for each column :");
    disp(norm(A,"inf",1));

    disp("Infinite norm of matrix values for each row :");
    disp(norm(A,"inf",2));

    
    disp("+----------------+");
    disp("| MATRIX PRODUCT |");
    disp("+----------------+");

    disp("Matrix-matrix product (C = A*B) :");
    disp(mtimes(ones<int>(3,4),eye<logical>(4,4)));
    
    disp("Blas form of matrix-matrix product (C = alpha*A*B + beta*C) :");
    int m=3, n=4;
    matrix<float> D = ones<float>(m,n);
    tgemm(M_PI,ones<int>(m,n),eye<logical>(n,n),2,D);
    disp(D);
    
    disp("+--------------------+");
    disp("| KROENECKER PRODUCT |");
    disp("+--------------------+");

    disp(kron(eye<int>(2,3),ones<float>(3,2)));
    
    
    disp("+--------------------+");
    disp("|    CROSS PRODUCT   |");
    disp("+--------------------+");

    disp(cross(eye(4,3),ones(4,3)));

    
    disp("+--------------------+");
    disp("|     DOT PRODUCT    |");
    disp("+--------------------+");

    disp(dot(eye(4,3),ones(4,3)));
    
    disp("+--------------------+");
    disp("|       UNIQUE       |");
    disp("+--------------------+");
    matrix<> B;
    matrix<std::size_t> Ia, Ib;
    A = eye(3);
    B = unique(A);
    std::tie(Ia,Ib) = argunique(eye(3));
    disp(eval(A(Ia)));
    disp(eval(B(Ib)));

    disp("+--------------------+");
    disp("|     INTERSECT      |");
    disp("+--------------------+");
    matrix<> C;
    A = eye(3);
    B = zeros(4);
    C = intersect(A,B);
    std::tie(Ia,Ib) = argintersect(A,B);
    disp(eval(A(Ia)));
    disp(eval(B(Ib)));
    
    disp("+--------------------+");
    disp("|      SETDIFF       |");
    disp("+--------------------+");
    A = eye(3);
    B = zeros(4);
    C = setdiff(A,B);
    Ia = argsetdiff(A,B);
    disp(eval(A(Ia)));
    
    disp("+--------------------+");
    disp("|       UNION        |");
    disp("+--------------------+");
    disp(union2(eye(3),zeros(4)));
    
    disp("+--------------------+");
    disp("|       GMRES        |");
    disp("+--------------------+");
    A = eye(100) + 1e-1*rand(100);
    B = rand(100,1);
    C = gmres(A,B);
    disp( norm(mtimes(A,C)-B,"inf") );
    C = gmres(A,B,1e-3,10);
    disp( norm(mtimes(A,C)-B,"inf") );
    matrix<> Am1 = gmres(A,eye(100));
    C = gmres(A,B,1e-3,10,Am1,C);
    disp( norm(mtimes(A,C)-B,"inf") );
    
    disp("+--------------------+");
    disp("|     CONVOLUTION    |");
    disp("+--------------------+");
    disp(conv(ones(3,4),eye(3),2));
    
    disp("+--------------------------+");
    disp("|     FOURIER TRANSFORM    |");
    disp("+--------------------------+");
    disp(dft(eye(1,3)));
    disp(real(idft(ones<std::complex<double>>(1,3))));
    
    
    disp("done !");
    return 0;
}
