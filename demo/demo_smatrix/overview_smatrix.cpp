/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : overview_smatrix.cpp                          |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Overview of sparse matrix operations          |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/smatrix.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|      INITIALIZE     |" << std::endl;
    std::cout << "+=====================+" << std::endl;
    
    // Documentation
    documentationFiles =
    {
        "/usr/local/include/castor/matrix.hpp",
        "/usr/local/include/castor/smatrix.hpp"
    };
    help("smatrix");
    
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|       INTERNAL      |" << std::endl;
    std::cout << "+=====================+" << std::endl;
    
    // Constructors
    smatrix<int> s;
    disp(s);
    s = 1;
    disp(s);
    smatrix<> As(3,4);
    disp(As);
    As = smatrix<>(3,4,5);
    disp(As);
    As = smatrix<>(2,3,{1,0,3,2,5,4},std::vector<double>({0.,0.,1.,2.,0.,3.}));
    disp(As);
    As = smatrix<>(2,3,{0,0,1,1},{0,1,0,1},std::vector<double>({1,0,0,2}));
    disp(As);
    As = M_PI*eye(3,4);
    disp(As);
    smatrix<float> Bs = As;
    disp(Bs);

    // Operators
    As = speye(2,3);
    disp(As);
    As += 1;
    disp(As);
    As += spones(2,3);
    disp(As);
    As -= 1;
    disp(As);
    As -= spones(2,3);
    disp(As);
    As *= 2;
    disp(As);
    As *= 2*spones(2,3);
    disp(As);
    As /= 2;
    disp(As);
    As /= 2*spones(2,3);
    disp(As);

    // As(l) and As(i,j) const
    smatrix<>const Fs(3,3,{1,1},{0,2},std::vector<double>({1,1}));
    disp(Fs);
    disp(Fs(0,0));
    disp(Fs(2,2));
    disp(Fs(1,0));
    disp(Fs(1,1));
    disp(Fs(1,2));
    disp(Fs);

    // As(l) and As(i,j) non const
    Bs = zeros(3,3);
    Bs(1,0) = 1;
    Bs(1,2) = 3;
    Bs(1,1) = 2;
    Bs(0,0) = 1;
    Bs(0,1) = 2;
    Bs(2,1) = 2;
    Bs(2,2) = 3;
    disp(Bs,2,std::cout,10);
    disp(Bs(0,2));
    disp(Bs(1,0));
    disp(Bs,2,std::cout,10);

    // As.ind(k) and As.val(k);
    Bs = speye(3,4);
    disp(Bs.ind(2));
    disp(Bs.val(2));

    // Tools
    int un = (int)s;
    disp(un);
    Bs = reshape(colon(1,12),3,4);
    Bs.reshape(4,3);
    disp(Bs);
    Bs.clear();
    disp(Bs);
    Bs = zeros(2,2);
    Bs(2) = 0; Bs(0) = 1; Bs(3) = 2; Bs(1) = 0;
    disp(Bs);
    check(Bs);
    disp(Bs);

    // Complex
    smatrix<std::complex<double>> Zs = ones(2,3) + M_1I*ones(2,3);
    Zs(1) = 0;
    disp(Zs);

        //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|       EXTERNAL      |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    // Builders
    disp(speye(3,4));
    disp(spones(3,4));
    disp(sprand(3,4));
    disp(spzeros(3,4));

    // Transformation
    disp(full(speye(3,4)));
    disp(sparse(eye(3,4)));
    smatrix<> Ms = reshape(colon(1,12),3,4);
    disp(full(reshape(Ms,4,3)));
    disp(full(transpose(Ms)));
    clear(Ms);
    disp(Ms);

    // Operator term-by-term
    disp(1+speye(2,3));
    disp(speye(2,3)+1);
    disp(speye(2,3)+ones(2,3));
    disp(ones(2,3)+speye(2,3));
    disp(speye(2,3)+spones(2,3));
    disp(-speye(2,3));
    disp(1-speye(2,3));
    disp(speye(2,3)-1);
    disp(speye(2,3)-spones(2,3));
    disp(2*speye(2,3));
    disp(speye(2,3)*2);
    disp(speye(2,3)*(2*spones(2,3)));
    disp(2/speye(2,3));
    disp(speye(2,3)/2);
    disp(speye(2,3)/(2*spones(2,3)));

    // Sparse Full product
    matrix<> A=rand(4,3), B=rand(3,4);
    A = vertcat(A,zeros(1,3));
    B = horzcat(zeros(3,1),B);
    matrix<> C=mtimes(A,B);
    disp(norm(mtimes(sparse(A),B)-C,"inf"));
    disp(norm(mtimes(A,sparse(B))-C,"inf"));
    disp(norm(mtimes(sparse(A),sparse(B))-C,"inf"));

    // Sparse to (bi-) linear index and values
    Ms = speye(3);
    disp(values(Ms));
    disp(index(Ms));
    matrix<std::size_t> I, J;
    matrix<> V;
    As = speye(3,4);
    As(2,3) = 1;
    std::tie(I,J,V) = find(As);
    disp(full(sparse(I,J,V)));
    
    // Gmres
    As = rand(5);
    B  = rand(5,2);
    matrix<> X = gmres(As,B,1e-6,10,speye(5),B);
    disp(norm(mtimes(As,X)-B,"inf"));
    
    // Perfo
    As = eye(1e2);
    Bs = rand(1e2);
    tic();
    smatrix<> Ds = mtimes(As,Bs);
    toc();
    disp(norm(full(Ds)-full(Bs),"inf"));
    
    disp("done !");
    return 0;
}
