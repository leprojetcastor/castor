/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : internal.cpp                                  |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Internal operators of matrix class            |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+======================+");
    disp("|  INTERNAL OPERATORS  |");
    disp("+======================+");
    
    matrix<> V = {0,1,2,3};
    matrix<float> A = eye(3,4);
    
    disp("operator+=");
    std::cout << "V =      : " << V << std::endl;
    V += 1;
    std::cout << "V += 1   : " << V << std::endl;
    V += {1,1,1,1};
    std::cout << "V += {}  : " << V << std::endl;
    V += std::vector<int>(4,1);
    std::cout << "V += v   : " << V << std::endl;
    V += matrix<double>(1,4,1);
    std::cout << "V += A   : " << V << std::endl;
    V += matrix<int>(1,4,1);
    std::cout << "V += A<> : " << V << std::endl;
    
    disp("operator-=");
    V = {0,1,2,3};
    std::cout << "V =      : " << V << std::endl;
    V -= 1;
    std::cout << "V -= 1   : " << V << std::endl;
    V -= matrix<int>(1,4,1);
    std::cout << "V -= A<> : " << V << std::endl;
    
    disp("operator*=");
    V = {0,1,2,3};
    std::cout << "V =      : " << V << std::endl;
    V *= 2;
    std::cout << "V *= 2   : " << V << std::endl;
    V *= matrix<int>(1,4,2);
    std::cout << "V *= A<> : " << V << std::endl;
    
    disp("operator/=");
    V = {0,1,2,3};
    std::cout << "V =      : " << V << std::endl;
    V /= 2;
    std::cout << "V /= 2   : " << V << std::endl;
    V /= matrix<int>(1,4,2);
    std::cout << "V /= A<> : " << V << std::endl;

    disp("operator()");
    V = {0,1,2,3};
    std::cout << "V         : " << V << std::endl;
    std::cout << "V(0)      : " << V(0) << std::endl;
    V(0) = -1;
    std::cout << "V(0)=-1   : " << V << std::endl;
    std::cout << "V(0,0)    : " << V(0,0) << std::endl;
    V(0,0) = -2;
    std::cout << "V(0,0)=-2 : " << V << std::endl;
    
    disp("Container of values stored in matrix A (std::vector<>) :");
    std::vector<float>const& val = A.val();
    for (std::size_t l=0; l<val.size(); ++l)
    {
        std::cout << val[l] << " ";
    }
    std::cout << std::endl;
       
    disp("Size of the matrix A :");
    disp(A.size());
    disp(A.size(1));
    disp(A.size(2));
    
    disp("Reshape matrix :");
    A.reshape(4,3);
    disp(A);

    disp("Resize matrix and fill with NAN :");
    A.resize(4,5);
    disp(A);
    
    disp("Clear matrix :");
    A.clear();
    disp(A);
    
    disp("done !");
    return 0;
}
