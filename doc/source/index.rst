Welcome to the castor library documentation
===========================================

The objective of the **castor** library is to propose **high-level semantics**, inspired by the Matlab language, allowing fast software prototyping in a low-level compiled language. It is nothing more than a **matrix management layer** using the tools of the standard **C++** library, in different storage formats (full, sparse and hierarchical). Indeed, the use of IDEs 1 such as Xcode, Visual studio, Eclipse, etc. allows today to execute compiled code (C, C++, fortran, etc.) with the **same flexibility as interpreted languages** (Matlab, Python, Julia, etc.). 

A **header-only** template library for matrix management has been developed based on the standard C++ library, notably the std::vector class. Many tools and algorithms are provided to simplify the development of scientific computing programs. Particular attention has been paid to semantics, for a simplicity of use "à la matlab", but written in C++. This high-level semantic/low-level language coupling makes it possible to gain efficiency in the prototyping phase, while ensuring performance for applications. In addition, direct access to data allows users to optimize the most critical parts of their code in native C++. Finally, **complete documentation** is available, as well as continuous integration unit tests. All of this makes it possible to meet the needs of teaching, academic issues and industrial applications at the same time. 

The **castor** library provides tools to : 

- create and manipulate dense, sparse and hierarchical matrices
- make linear algebra computations based on optimized BLAS library
- make graphical representations based on VTK library

These tools are used by applicative projects : 

- finite and boundary element method using Galerkin approximation
- analytical solutions for scattering problems 

The source files of the library is available here : `<https://gitlab.labos.polytechnique.fr/leprojetcastor/castor>`_.

As the semantics offered by **castor** library being voluntarily close to the Matlab environment, there are functions signature and their documentation inspired by it. You can refer to https://fr.mathworks.com/help/matlab/index.html.

Licensing
---------

The **castor** library is provided in open source under LGPL 3.0.

This program is free software, distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,  redistribute and/or modify it under the terms of the GNU Lesser General Public License, as published by the Free Software Foundation (version 3 or later,  http://www.gnu.org/licenses).


Gallery
-------
.. image:: img/head.png
   :width: 600
   :align: center

Boundary element computation using le projet castor. Simulation of acoustic scattering by human head, excited by plane wave at 8Khz (20.000 degrees of freedom). Hierarchical solver is used on a standard laptop.

.. toctree::
   :caption: Installation
   :maxdepth: 1
   :hidden:

   installation

.. _label-user-guide:

.. toctree::
   :caption: User guide
   :maxdepth: 1
   :hidden:

   getting_started
   basics
   advanced
   linalg
   graphics
   sparse_matrix 
   kissfft

.. toctree::
   :caption: Dense matrix
   :maxdepth: 1
   :hidden:

   class_matrix
   class_view_cview
   algorithm
   builders
   geometry
   io
   math
   dimensions
   manipulations
   operators
   tools 
   transforms 

.. toctree::
   :caption: Linear algebra
   :maxdepth: 1
   :hidden:

   factorization
   linear_solver
   singular_eig_values
   lowlevel_linalg_func

.. toctree::
   :caption: Graphical rendering
   :maxdepth: 1
   :hidden:

   class_figure
   basic_plot
   graphical_io
   graphical_tools
   mesh_management
   mesh_plot

.. toctree::
   :caption: Sparse matrix
   :maxdepth: 1
   :hidden:

   class_smatrix
   api_smatrix 

.. toctree::
   :caption: Hierarchical matrix
   :maxdepth: 1
   :hidden:

   class_hmatrix
   class_bintree
   api_hmatrix 

.. toctree::
   :caption: Applications
   :maxdepth: 1
   :hidden:

   fembem
   analyticalscattering

.. toctree::
   :caption: Contacts
   :maxdepth: 1
   :hidden:
 
   developpers

