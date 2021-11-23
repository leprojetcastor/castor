/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : demo_smatrix.cpp                          |
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
    smatrix<> As = {1,2,3};
    disp(As);
    As = {{1,2,3},{4,5,6}};
    disp(As);
    As = smatrix<>(3,4);
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
    check(Zs);
    disp(Zs);

    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|       EXTERNAL      |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    // Builders (m,n)
    disp(speye(3,4));
    disp(spones(3,4));
    disp(sprand(3,4));
    disp(spzeros(3,4));
    
    // Builders {m,n}
    disp(speye({3,4}));
    disp(spones({3,4}));
    disp(sprand({3,4}));
    disp(spzeros({3,4}));
    
    // Dimensions
    disp(size(speye(3,4)));
    disp(size(speye(3,4),1));
    disp(size(speye(3,4),2));
    disp(numel(speye(3,4)));
    disp(length(speye(3,4)));
    disp(nnz(speye(3,4)));

    // Transformation
    disp(full(speye(3,4)));
    disp(sparse(eye(3,4)));
    smatrix<> Ms = reshape(colon(1,12),3,4);
    disp(full(reshape(Ms,4,3)));
    disp(full(transpose(Ms)));
    clear(Ms);
    disp(Ms);

    // Operator+
    disp(1+speye(2,3));
    disp(speye(2,3)+1);
    disp(speye(2,3)+ones(2,3));
    disp(ones(2,3)+speye(2,3));
    disp(speye(2,3)+spones(2,3));

    // Operator-
    disp(-speye(2,3));
    disp(1-speye(2,3));
    disp(speye(2,3)-1);
    disp(speye(2,3)-ones(2,3));
    disp(ones(2,3)-speye(2,3));
    disp(speye(2,3)-spones(2,3));

    // Operator*
    disp(2*speye(2,3));
    disp(speye(2,3)*2);
    disp(speye(2,3)*(2*ones(2,3)));
    disp((2*ones(2,3))*speye(2,3));
    disp(speye(2,3)*(2*speye(2,3)));

    // Operator/
    disp(2/speye(2,3));
    disp(speye(2,3)/2);
    disp(speye(2,3)/(2*ones(2,3)));
    disp((2*ones(2,3))/speye(2,3));
    disp(speye(2,3)/(2*speye(2,3)));

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
    matrix<> V;
    As = speye(3,4);
    As(2,3) = 1;
    matrix<std::size_t> I, J;
    std::tie(I,J,V) = find(As);
    disp(full(sparse(I,J,V)));
    
    // Gmres
    As = speye(100) + 1e-3*rand(100);
    B  = rand(100,4);
    matrix<> X = gmres(As,B);
    disp(norm(mtimes(As,X)-B,"inf"));
    smatrix<> Asm1 = gmres(As,eye(100));
    X = gmres(As,B,1e-3,10,Asm1,X);
    disp(norm(mtimes(As,X)-B,"inf"));
    
    // Linsolve
    X = linsolve(As,B);
    disp(norm(mtimes(As,X)-B,"inf"));
    
    // Spdiags
    matrix<> e = reshape(colon(1,9),3,3);
    disp(full(spdiags(e,{-2,0,2},3,6)));
    disp(full(spdiags(e,{-2,0,2},3,3)));
    disp(full(spdiags(e,{-2,0,2},6,3)));

    //===============================================================
    std::cout << "+=================+" << std::endl;
    std::cout << "|       VIEW      |" << std::endl;
    std::cout << "+=================+" << std::endl;
    
    matrix<std::size_t> L;
    I = {0,1,0,1};
    J = {1,2,2,1,0};
    L = {{5,2,0},{5,4,5}};
    
    As = matrix<>({{0,1,2},{3,4,5}});
    disp(full(As));

    Bs = get(As,L);
    disp(full(Bs));

    Bs = get(As,I,J);
    disp(full(Bs));
    
    As = spones(2,3);
    L  = all(As);
    set(As,L,2*speye(1,6));
    disp(full(As));

    As = spones(2,3);
    I  = row(As);
    J  = col(As);
    set(As,I,J,2*speye(2,3));
    disp(full(As));

    As = spones(2,3);
    Bs = 2*speye(2,3);
    disp(eval(Bs(L)));
    disp(eval(Bs(I,J)));

    As    = spones(2,3);
    As(L) = eval(Bs(L));
    disp(As);

    As      = spones(2,3);
    As(I,J) = eval(Bs(I,J));
    disp(As);
    
    smatrix<>const Gs = 2*speye(2,3);
    disp(eval(Gs(L)));
    disp(eval(Gs(I,J)));
    
    As = ones(4,4);
    As({3,1,2},{3,1,2}) = smatrix<>(diag<double>({1,2,3}));
    disp(full(As));
    disp(full(eval(As({3,1,2},{3,1,2}))));
    
    //===============================================================
    std::cout << "+=================+" << std::endl;
    std::cout << "|      PERFO      |" << std::endl;
    std::cout << "+=================+" << std::endl;
    
    As = eye(1e2);
    Bs = rand(1e2);
    tic();
    smatrix<> Ds = mtimes(As,Bs);
    toc();
    disp(norm(full(Ds)-full(Bs),"inf"));
    
    
    disp("done !");
    return 0;
}
