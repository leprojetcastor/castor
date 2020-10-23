/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : overview_hmatrix.cpp                          |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Overview of all fonctionalities of hmatrix    |
 |  `---'  |                framework                                     |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/linalg.hpp"
#include "castor/hmatrix.hpp"

using namespace castor;

template<typename T, typename S>
matrix<T> Gxy(matrix<S>const& X, matrix<S>const& Y, S wavenumber=0, matrix<std::size_t> Ix={}, matrix<std::size_t> Iy={});

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
        "/usr/local/include/castor/hmatrix.hpp"
    };
    help("");
    
    // Parameters
    std::size_t Nx  = 999;
    std::size_t Ny  = 999;
    double      tol = 1e-3;
    
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|     BINARY TREE     |" << std::endl;
    std::cout << "+=====================+" << std::endl;
    
    // Sphere X
    matrix<double> x, y, z, X;
    std::tie(x,y,z) = sphere2(Nx);
    X = horzcat(reshape(x,Nx,1),reshape(y,Nx,1));
    X = horzcat(X,reshape(z,Nx,1));
    
    // Sphere Y
    matrix<double> Y;
    std::tie(x,y,z) = sphere2(Ny);
    Y = horzcat(reshape(x,Ny,1),reshape(y,Ny,1));
    Y = horzcat(Y,reshape(z,Ny,1));
    //    Y = 1 + Y;
    
    // Tree
    binarytree<double> Xtree(X);
    binarytree<double> Ytree(Y);
    
    //===============================================================
    std::cout << "+==================+" << std::endl;
    std::cout << "|  TOTAL PIVOTING  |" << std::endl;
    std::cout << "+==================+" << std::endl;
    
    // Full matrix
    tic();
    matrix<double> M = Gxy<double>(X,Y);
    toc();
    
    // H-Matrix with total pivoting
    tic();
    hmatrix<double> Mh(X,Y,tol,M);
    toc();
    disp(Mh);
    disp( norm(M-full(Mh),"inf")/norm(M,"inf") );
    
    //===============================================================
    std::cout << "+====================+" << std::endl;
    std::cout << "|  PARTIAL PIVOTING  |" << std::endl;
    std::cout << "+====================+" << std::endl;

    // Define kernel function as lambda function
    auto fct = [&X,&Y](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return Gxy<double>(X,Y,0.,Ix,Iy);
    };

    // H-Matrix with partial pivoting
    tic();
    hmatrix<double> Mh2(X,Y,tol,fct);
    toc();
    disp(Mh2);
    disp( norm(M-full(Mh2),"inf")/norm(M,"inf") );

    //===============================================================
    std::cout << "+=================+" << std::endl;
    std::cout << "|     ALGEBRA     |" << std::endl;
    std::cout << "+=================+" << std::endl;
    
    // Invertible matrix
    M  = Gxy<double>(X,Y);
    M += sqrt(Nx)*eye(Nx,Ny);
    Mh = hmatrix<double>(X,Y,tol,M);
    disp( norm(M-full(Mh),"inf")/norm(M,"inf") );

    // FULL: Left product
    matrix<double> V   = rand(2,Nx);
    matrix<double> ref = mtimes(V,M);
    matrix<double> sol = mtimes(V,Mh);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // FULL: Right product
    V   = rand(Ny,2);
    ref = mtimes(M,V);
    sol = mtimes(Mh,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // SCALAR: Add
    ref = mtimes(1.+M+1.,V);
    sol = mtimes(1.+Mh+1.,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
  
    // SCALAR: Product
    ref = mtimes(M_PI*M*M_PI,V);
    sol = mtimes(M_PI*Mh*M_PI,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
    
    // SCALAR: Uminus
    ref = mtimes(-M,V);
    sol = mtimes(-Mh,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
    
    // SCALAR: Substract
    ref = mtimes(1.-M-1.,V);
    sol = mtimes(1.-Mh-1.,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // SCALAR: Right divide
    ref = mtimes(M/M_PI,V);
    sol = mtimes(Mh/M_PI,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // FULL: Add low-rank
    ref = mtimes(2+3.*M,V);
    hmatrix<double> tmp = Mh;
    tgeabm(2.,ones(Nx,1),ones(1,Ny),3.,tmp);
    sol = mtimes(tmp,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX : Transpose
    ref = mtimes(transpose(V),transpose(M));
    sol = mtimes(transpose(V),transpose(Mh));
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX: Add h-matrix
    ref = mtimes(2.*M+1,V);
    tmp = hmatrix<double>(Nx,Ny,tol,1);
    sol = mtimes(Mh+Mh+tmp,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
    
    // HMATRIX: TGEMM
    matrix<> tmp2 = M;
    tgemm(1,M,M,1,tmp2);
    ref = mtimes(tmp2,V);
    tmp = Mh;
    tgemm(1.,Mh,Mh,1.,tmp);
    sol = mtimes(tmp,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX: Multiply h-matrix
    ref = mtimes(mtimes(M,M),V);
    sol = mtimes(mtimes(Mh,Mh),V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
    
    // HMATRIX: Inverse
    hmatrix<double> Mhm1;
    ref = mtimes(inv(M),V);
    Mhm1 = inv(Mh);
    sol = mtimes(Mhm1,V);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
    disp(Mhm1);

    // HMATRIX: LU
    hmatrix<double> Lh, Uh;
    std::tie(Lh,Uh) = lu(Mh);
    ref = mtimes(M,V);
    sol = mtimes(Lh,mtimes(Uh,V));
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );
    disp(Lh);
    disp(Uh);
    
    // HMATRIX: LLOWSOLVE
    ref = mtimes(inv(Lh),V);
    sol = V;
    Lh.hllowsolve(sol);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX: RUPSOLVE
    ref = mtimes(transpose(V),inv(Uh));
    sol = transpose(V);
    Uh.hrupsolve(sol);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX: LUPSOLVE
    ref = mtimes(inv(Uh),V);
    sol = V;
    Uh.hlupsolve(sol);
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX: LINSOLVE
    tic();
    ref = linsolve(M,V);
    toc();
    tic();
    sol = linsolve(Mh,V);
    toc();
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    // HMATRIX: GMRES
    tic();
    ref = gmres(M,V,tol,100);
    toc();
    tic();
    sol = gmres(Mh,V,tol,100);
    toc();
    disp( norm(ref-sol,"inf")/norm(ref,"inf") );

    
    
    disp("done !");
    return 0;
}


// Green kernel function
template<typename T, typename S>
matrix<T> Gxy(matrix<S>const& X, matrix<S>const& Y, S wavenumber, matrix<std::size_t> Ix, matrix<std::size_t> Iy)
{
    // Initialize dof data
    if (numel(Ix)==0) {Ix = range(0,size(X,1));}
    if (numel(Iy)==0) {Iy = range(0,size(Y,1));}
    if (max(Ix)>size(X,1) || max(Iy)>size(Y,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Indices are not in data bound.");
    }
    matrix<T> M(numel(Ix),numel(Iy));
    
    // Initialize local data
    std::size_t ix, iy;
    S rxy;
    T imk=std::sqrt(T(-1))*wavenumber;
    
    // Build
    for (std::size_t i=0; i<numel(Ix); ++i)
    {
        ix = Ix(i);
        for (std::size_t j=0; j<numel(Iy); ++j)
        {
            iy  = Iy(j);
            rxy = std::sqrt((X(ix,0)-Y(iy,0))*(X(ix,0)-Y(iy,0)) +
                            (X(ix,1)-Y(iy,1))*(X(ix,1)-Y(iy,1)) +
                            (X(ix,2)-Y(iy,2))*(X(ix,2)-Y(iy,2)) );
            if (wavenumber==0)
            {
                M(i,j) = (rxy>1e-12)? (1./rxy) : (0);
            }
            else
            {
                M(i,j) = (rxy>1e-12)? (std::exp(imk*rxy)/rxy) : (imk);
            }
        }
    }
    return M;
}
