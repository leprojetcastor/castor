---
title: 'Castor: A C++ library to code "à la matlab"'
tags:
  - C++
  - Scientific computing
  - Fast prototyping
  - FEM/BEM simulation
authors:
  - name: Matthieu Aussal^[matthieu.aussal@polytechnique.edu] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-2812-7578
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Marc Bakry^[marc.bakry@polytechnique.edu] # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Laurent Series^[laurent.series@polytechnique.edu]
    affiliation: 2
affiliations:
 - name: Ecole Polytechnique (CMAP), INRIA, Institut Polytechnique Paris, Route de Saclay 91128, Palaiseau, France
   index: 1
 - name: Ecole Polytechnique (CMAP), Institut Polytechnique Paris, Route de Saclay 91128, Palaiseau, France
   index: 2
date: 22 October 2021
bibliography: paper.bib
---

# Summary

The objective of the *Castor* framework is to propose high-level semantics, inspired by the Matlab language, allowing fast software prototyping in a low-level compiled language. It is nothing more than a matrix management layer using the tools of the standard C++ library (C++14 and later), in different storage formats (full, sparse and hierarchical). Indeed, the use of IDEs such as Xcode, Visual studio, Eclipse, etc. allows today to execute compiled code (C, C++, fortran, etc.) with a flexibility close to the one of interpreted languages (Matlab, Python, Julia, etc.). The *Castor* framework is provided as an open source software under the LGPL 3.0.

# Statement of need

Matlab is a software used worldwide in numerical prototyping, due to its particularly user-friendly semantics and its certified toolboxes. However, some usecases do not allow codes in Matlab format, for example multi-platform portability issues or proprieraty licensing. To start meeting these needs, a header-only template library for matrix management has been developed, based on the standard C++ library, by encapsulating the `std::vector` class. Many tools and algorithms are provided to simplify the development of scientific computing programs:
 
 - dense, sparse and hierarchical matrices manipulations,
 - linear algebra computations (optimized BLAS library),
 - graphical representations (VTK library).

Compared to standard C++, this high-level semantic/low-level language coupling makes it possible to gain efficiency in the prototyping phase, while ensuring performance for applications. In addition, direct access to data structures allows users to optimize the most critical parts of their code. Finally, a complete documentation is available, as well as continuous integration unit tests. All of this makes it possible to meet the needs of teaching, academic issues and industrial applications at the same time.


# Dense Matrix  

The dense matrix part of the *Castor* framework implements its own templatized class `matrix<T>` in *matrix.hpp* , where `T` can be fundamental types of C++ as well as `std::complex`. This class is built over a `std::vector<T>` which holds the values [@cpp:14]. Note that the element of a matrix is stored in row-major order and that a vector is considered as a $1\times n$ or $n\times 1$ matrix.

The class `matrix<T>` provides many useful functions and operators such as:

- builders which can be used to initialize all coefficients (`zeros`, `ones`, `eye`, etc.),
- standard algorithms over data stored in matrices (`norm`, `max`, `sort`, `argsort`, etc.),
- mathematical functions which can be applied element-wise (`cos`, `sqrt`, `conj`, etc.),
- matrix manipulations like concatenate matrices in all dimensions, find the non-zero elements or transpose them, reshape size, etc.,
- standard C++ operators which have been overloaded and work element-wise (`+`,`*`,`!`,`&&`,etc.),
- values accessors and matrix views with linear and bi-linear indexing,
- elements of linear algebra, such as the matrix product (`mtimes` or `tgemm`) and linear system resolution (multi-right-hand-side `gmres`),
- many other tools to display elements (`<<`, `disp`), save and load elements from file (ASCII or binary), etc.

The API provides more than a hundred functions and is designed such that it should feel like using Matlab. For advanced users, direct access to the data stored in the `std::vector<T>` enables all or part of an algorithm to be optimized in native C++. 

This example displays the sum of two matrices with implicit cast :

```c++
#include <iostream>
#include "castor/matrix.hpp"
using namespace castor;
int main(int argc, char* argv[])
{
    matrix<float> A = {{ 1.0,  2.0,  3.0,  4.0},
                       { 5.0,  6.0,  7.0,  8.0},
                       { 9.0, 10.0, 11.0, 12.0}};
    matrix<double> B = eye(3,4);
    auto C = A + B;
    disp(C);
    return 0;
}
```
```
Matrix 3x4 of type 'd' (96 B):
    2.00000      2.00000      3.00000      4.00000  
    5.00000      7.00000      7.00000      8.00000  
    9.00000     10.00000     12.00000     12.00000  
```

# Linear Algebra

The linear algebra part of the framework, implemented in *linalg.hpp*, provides a set of useful functions to perform linear algebra computations by linking to optimized implementations of the BLAS and LAPACK standards [@blaslapack:00] (OpenBLAS, oneAPI MKL, etc.).

The BLAS part is a straightforward overlay of the C-BLAS type III API, which is compatible with row-major ordering.  This is achieved by a template specialization of the `tgemm` function, which allows optimized implementations to take control of the computation using `sgemm`, `dgemm`, `cgemm` and `zgemm`. Thanks to this interface, naive implementations proposed in *matrix.hpp* for dense matrix-products `mtimes` and `tgemm` may be improved in term of performance, especially for large matrices. 

The LAPACK part is a direct overlay over the Fortran LAPACK API, which uses a column ordering storage convention. This interface brings new high-level functionalities, such as a linear solver (`linsolve`), matrix inversion (`inv`, `pinv`), factorizations (`qr`, `lu`), the search for eigen or singular values decompositions (`eig` ,`svd`), aca compression (`aca`), etc. It uses templatized low-level functions following the naming convention close to the LAPACK one (like `tgesdd`, `tgeqrf`, etc.).

This example displays the product of $A$ and $A^{-1}$ :

```c++
#include <iostream>
#include "castor/matrix.hpp"
#include "castor/linsolve.hpp"
using namespace castor;
int main (int argc, char* argv[])
{
    matrix<> A = rand(4);
    matrix<> Am1 = inv(A); 
    disp(mtimes(A,Am1));
    return 0;
}
```
```
Matrix 4x4 of type 'd' (128 B):
 1.0000e+00   1.0408e-16  -2.7756e-17  -5.5511e-17  
          0   1.0000e+00  -5.5511e-17   1.1102e-16  
          0  -2.2204e-16   1.0000e+00  -1.1102e-16  
-2.7756e-17            0            0   1.0000e+00  
```

# 2D/3D Visualization

The graphic rendering part, provided by *graphics.hpp*, features 2D/3D customizable plotting and basic mesh generation. It is based on the well-known VTK library [@vtk:2000]. Here again, the approach tries to get as close as possible to Matlab semantics.

First, the user creates a `figure`, which is a dynamic container of data to display. The `figure` class is composed of a `vtkContextView` class, providing a view with a default interactor style, renderer, etc. Then, graphic representations can be added to the figure, using functions like `plot`, `imagesc`, `plot3`, `mesh`, etc. Options are available to customize the display of the results, such as the plotting style, legend, colorbar and others basic stuff. Finally, the `drawnow` function must be called to display all defined figures. The latters are displayed and manipulated in independent windows.

In addition, graphics exports are available in different compression formats (`png`,` jpg`, `tiff`, etc.), as well as video rendering (`ogg`). 

This example shows a basic 2D plotting of a sine function:

```c++
#include "castor/matrix.hpp"
#include "castor/graphics.hpp" 
using namespace castor;
int main (int argc, char* argv[])
{
    matrix<> X = linspace(0,10,100);
    figure fig;
    plot(fig,X,sin(X),{"r-+"},{"sin(x)"});
    plot(fig,X,cos(X),{"bx"},{"cos(x)"});
    drawnow(fig);
    return 0;
}
```
![Caption for example figure.\label{fig:sin}](plot2d.png)

# Sparse matrices

Some matrices have sparse structures, with many (many) zeros that do not need to be stored [@sparse:2011]. There are adapted storage formats for this type of structure (LIL, COO, CSR, etc.), the most natural being to store the indices of rows and columns for each non-zero value, as a list of triplet $\{i,j,v\}$. For the *Castor* framework, a dedicated template class to this kind of matrix has been developed (see *smatrix.hpp*). The storage format is based on a row major sorted linear indexing. Only non-zero values and their sorted linear indices are stored in a list of pairs $\{v,l\}$:  for a $m\times n$ matrix, the following bijection is used to switch with the common bilinear indexation:

 - $\{i,j\} \rightarrow l = i \cdot n + j$, 
 - $l \rightarrow \{i=\frac{l}{n}; j= i\textrm{ mod }n\}$. 

Accessors to all the elements are provided so that sparse matrices can be manipulated in a similar way as the dense matrices. This operation is performed by dichotomy with a convergence in $log_2 (\text{nnz})$, where $\text{nnz}$ is the number of non-zero elements. Just like dense matrices, numerical values are stored in a templatized `std::vector<T>`. For convenience, we provide classical builders (`sparse`, `speye`, `spdiags`, etc.), standard C++ operators overloading, views, display functions (`disp`, `spy`) and some linear algebra tools (`transpose`, `mtimes`, `gmres`, etc.). 

This example displays the sum of two sparse matrices, with implicit cast and sparse to dense conversion :

```c++
#include <iostream>
#include "castor/smatrix.hpp"
using namespace castor;
int main (int argc, char* argv[])
{
    smatrix<float> As = {{0.0,  0.0,  0.0},
                         {5.0,  0.0,  7.0}};
    As(0,1) = 2.0;
    smatrix<double> Bs = speye(2,3);
    disp(As);
    disp(As(0,1)); // bilinear accessor
    disp(As(4));   // linear accessor
    disp(Bs);
    disp(full(As+Bs));
    return 0;
}
```
```
Sparse matrix 2x3 of type 'f' with 3 elements (12 B):
(0,1)  2
(1,0)  5
(1,2)  7
2
0
Sparse matrix 2x3 of type 'd' with 2 elements (16 B):
(0,0)  1
(1,1)  1
Matrix 2x3 of type 'd' (48 B):
    1.00000      2.00000            0  
    5.00000      1.00000      7.00000  
```

# Hierarchical matrices

To widen the field of applications, the $\mathcal H$-matrix format, so-called hierachical matrices [@hackbush:1999], have been added in *hmatrix.hpp*. They are specially designed for matrices with localized rank defaults. It allows a fully-populated matrix to be assembled and stored in a lighter format by compressing some parts of the original dense matrix using a low-rank representation [@rjasanow:2002]. They are constructed by binary tree subdivisions in a recursive way, with a parallel assembly of the compressed and full leaves (using the OpenMP standard). This format features a complete algebra, from elementary operations to matrix inversion. An example is given in the application section that follows.

# Application with a FEM/BEM simulation
As an application example, an acoustical scattering simulation was carried out using a boundary element method (BEM) tool, implemented with the *Castor* framework (see the *fembem* package [@fembem:21]). We consider a smooth $n$-oriented surface $\Gamma$ of some object $\Omega$, illuminated by an incident plane wave $u_i$ with wave-number $k$. The scattered field $u$ satisfies the Helmholtz equation in $\Omega$, Neumann boundary conditions (*sound-hard*) and the Sommerfeld radiation condition: 

- $-(\Delta u + k^2 u) = 0$, 
- $-\partial_n u_i = 0$, 
- $\lim\limits_{r \to + \infty} r\textrm{ }(\partial_r u - i k u) = 0$. 
  
The scattered field $u$ satisfies the integral representation (Neumann interior extension, see [@terasse:2007]):

$$u(\textbf{x}) = - \int_\Gamma \partial_{n_y} G(\textbf{x},\textbf{y})\mu(\textbf{y}) d_y  \quad \forall  \textbf{x} \in \mathbb{R}^3 \setminus \overline \Omega,$$

for some density $\mu$, with the Green kernel $G(\textbf{x},\textbf{y}) = \displaystyle\frac{e^{i k |x - y|}}{4 \pi |x - y| }$. Using the boundary conditions we obtain :

$$ - H\mu(\textbf{x})  = - \partial_n u_i(\textbf{x}) \quad \forall \textbf{x} \in \Gamma,$$

where the hypersingular operator $H$ is defined by:

$$H \mu(\textbf{x}) = \int_\Gamma \partial_{n_x} \partial_{n_y} G(\textbf{x},\textbf{y})\mu(\textbf{y}) d_y.$$

The operator $H$ is assembled using a $P_1$ finite element discretization on a triangular mesh of the surface $\Gamma$, stored using dense matrices (*matrix.hpp*) or hierarchical matrices (*hmatrix.hpp*).

![Caption for example figure.\label{fig:head}](head.png)

Finaly, using all the tools provided by Castor to write and solve these equations, we are able to efficiently compute the acoustic diffraction of a harmonic plane wave at 8kHz, on a human head mesh [@symare:2013]). The simulation result (\autoref{fig:head}) highlights the role of the auditory pavilion as a resonator, modifying the timbre of a sound source to allow a listener's brain to precisely locate its direction. 

```c++
#include <castor/matrix.hpp>
#include <castor/smatrix.hpp>
#include <castor/hmatrix.hpp>
#include <castor/linalg.hpp>
#include <castor/graphics.hpp>
#include "fem.hpp"
#include "bem.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    // Load meshes
    matrix<double> Svtx;
    matrix<size_t> Stri;
    std::tie(Stri,Svtx) = triread("./","Head03_04kHz.ply");
    
    // Graphical representation
    figure fig;
    trimesh(fig,Stri,Svtx);
    
    // Parameters
    matrix<double> U = {0,0,-1};
    double f = 2000;
    double k = 2*M_PI*f/340;
    float tol = 1e-3;
    
    // FEM and mass matrix (sparse storage)
    tic();
    femdata<double> v(Stri,Svtx,lagrangeP1,3);
    femdata<double> u(Stri,Svtx,lagrangeP1,3);
    auto Id = mass<std::complex<double>>(v);
    toc();

    // Left hand side : [-H]
    tic();
    auto LHSfct = [&v,&u,&k](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return -hypersingular<std::complex<double>>(v,u,k,Ix,Iy);
    };
    hmatrix<std::complex<double>> LHS(v.dof(),u.dof(),tol,LHSfct);
    toc();
    disp(LHS);

    // Right hand side : - \int_sigma_x v(x) PW(x,u) dx
    auto B = - rightHandSide<std::complex<double>>(v,dnPWsource,U,k);

    // Solve -H = -dnP0
    hmatrix<std::complex<double>> Lh,Uh;
    tic();
    std::tie(Lh,Uh) = lu(LHS,1e-1);
    toc();
    disp(Lh);
    disp(Uh);
    auto mu = gmres(LHS,B,tol,100,Lh,Uh);
    
    // Boundary radiation
    tic();
    auto Dbndfct = [&v,&u,&k,&Id](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return 0.5*eval(Id(Ix,Iy)) + doubleLayer<std::complex<double>>(v,u,k,Ix,Iy);
    };
    hmatrix<std::complex<double>> Dbnd(v.dof(),u.dof(),tol,Dbndfct);
    matrix<std::complex<double>> Pbnd = mtimes(Dbnd,mu);
    toc();
    Pbnd = - gmres(Id,Pbnd,tol,100) + planeWave(v.dof(),U,k);
    
    // Graphical representation
    figure fig2;
    trimesh(fig2,Stri,Svtx,abs(Pbnd));
    
    // Export in .vtk file
    triwrite("./","head.vtk",Stri,Svtx,real(Pbnd));
    
    // Plot
    drawnow(fig);
    disp("done !");
    return 0;
}
```
  
# Acknowledgements

We thank Houssem Haddar for his precious help.

# References
