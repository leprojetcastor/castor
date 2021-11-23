/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : perf_hmat.cpp                              |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 14.02.2021                                    |
 |  / 0 \  |   LAST MODIF : 14.02.2021                                    |
 | ( === ) |   SYNOPSIS   : H-matrix performances and precision test      |
 |  `---'  |                for particular kernel                         |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/linalg.hpp"
#include "castor/smatrix.hpp"
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
    help("hmatrix");
    
    // Parameters
    std::size_t Nx  = 1e3;
    std::size_t Ny  = Nx;  // NUMBER OF EMMITERS
    std::size_t Nr  = 100;  // NUMBER OF EMMITERS
    double      tol = 1e-3; // ACCURACY
    double      f   = 100;  // FREQUENCY (Hz)
    double      c   = 343;  // SOUND CELERITY (m.s-1)
    
    // Wave numbers
    double k = 2*M_PI*f/c;
    
    // Sphere X
    disp("---> DISCRETIZATION X <---");
    matrix<double> x, y, z, X;
    std::tie(x,y,z) = sphere2(Nx);
    X = horzcat(reshape(x,Nx,1),reshape(y,Nx,1));
    X = horzcat(X,reshape(z,Nx,1));
    disp(bbox(X));
    
    // Sphere Y
    disp("---> DISCRETIZATION Y <---");
    matrix<double> Y;
    std::tie(x,y,z) = sphere2(Ny);
    Y = horzcat(reshape(x,Ny,1),reshape(y,Ny,1));
    Y = horzcat(Y,reshape(z,Ny,1));
    //        Y = -5 + 0.5*Y;
    disp(bbox(Y));
    disp(bbox(X).isfar(bbox(Y)));
    
    // Full computations
    if (Nx<=1e3)
    {
        // Full
        disp("---> FULL <---");
        tic();
        matrix<std::complex<double>> M = Gxy<std::complex<double>>(X,Y,k);
        toc();
        
        // H-Matrix from full
        disp("---> FULL TO H-MATRIX <---");
        tic();
        hmatrix<std::complex<double>> Mh(X,Y,tol,M);
        toc();
        disp(Mh);
        disp(norm(M-full(Mh),"inf")/norm(M,"inf"));
        
        // H-Matrix from sparse
        disp("---> SPARSE TO H-MATRIX <---");
        tic();
        Mh = hmatrix<std::complex<double>>(X,Y,tol,sparse(M));
        toc();
        disp(Mh);
        disp(norm(M-full(Mh),"inf")/norm(M,"inf"));
    }
    
    // FUll reference
    disp("---> FULL MATRIX REFERENCE <---");
    matrix<std::complex<double>>   V = rand(Ny,2) + M_1I*rand(Ny,2);
    matrix<std::size_t>            I = round(linspace(0,Nx-1,Nr));
    tic();
    matrix<std::complex<double>>  Mr = Gxy<std::complex<double>>(X,Y,k,I,range(0,Ny));
    matrix<std::complex<double>> MrV = mtimes(Mr,V);
    toc();
    
    // H-Matrix from lambda
    disp("---> H-MATRIX <---");
    auto fct = [&X,&Y,&k](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return Gxy<std::complex<double>>(X,Y,k,Ix,Iy);
    };
    tic();
    hmatrix<std::complex<double>> Mh(X,Y,tol,fct);
    toc();
    disp(Mh);
    matrix<std::complex<double>> MhV = mtimes(Mh,V);
    disp(norm(MrV-eval(MhV(I,col(V))),"inf")/norm(MrV,"inf"));
    
    // H-Matrix LU
    disp("---> FACTO LU <---");
    hmatrix<std::complex<double>> Lh, Uh;
    tic();
    std::tie(Lh,Uh) = lu(Mh);
    toc();
    disp(Lh);
    disp(Uh);
    matrix<std::complex<double>> LhUhV = mtimes(Lh,mtimes(Uh,V));
    disp(norm(MrV-eval(LhUhV(I,col(V))),"inf")/norm(MrV,"inf"));

    // GMRES precond LU
    disp("---> GMRES + PRECOND LU <---");
    tic();
    std::tie(Lh,Uh) = lu(Mh,1e-1);
    toc();
    disp(Lh);
    disp(Uh);
    LhUhV = mtimes(Lh,mtimes(Uh,V));
    disp(norm(MrV-eval(LhUhV(I,col(V))),"inf")/norm(MrV,"inf"));
    gmres(Mh,V,tol,100,Lh,Uh);
    
    // Exit
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
            M(i,j) = (rxy>1e-12)? (std::exp(imk*rxy)/rxy) : (1.);
        }
    }
    return M;
}



