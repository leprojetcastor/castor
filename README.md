The objective of the castor library is to propose high-level semantics, inspired by the Matlab language, allowing fast software prototyping in a low-level compiled language. It is nothing more than a matrix management layer using the tools of the standard C++ library, in different storage formats (full, sparse and hierarchical). Indeed, the use of IDEs 1 such as Xcode, Visual studio, Eclipse, etc. allows today to execute compiled code (C, C++, fortran, etc.) with the same flexibility as interpreted languages (Matlab, Python, Julia, etc.).

A header-only template library for matrix management has been developed based on the standard C++ library, notably the std::vector class. Many tools and algorithms are provided to simplify the development of scientific computing programs. Particular attention has been paid to semantics, for a simplicity of use “à la matlab”, but written in C++. This high-level semantic/low-level language coupling makes it possible to gain efficiency in the prototyping phase, while ensuring performance for applications. In addition, direct access to data allows users to optimize the most critical parts of their code in native C++. Finally, complete documentation is available, as well as continuous integration unit tests. All of this makes it possible to meet the needs of teaching, academic issues and industrial applications at the same time.

The castor library provides tools to :

    create and manipulate dense, sparse and hierarchical matrices

    make linear algebra computations based on optimized BLAS library

    make graphical representations based on VTK library

These tools are used by applicative projects :

    boundary element method using Galerkin approximation

    analytical solutions for scattering problems
