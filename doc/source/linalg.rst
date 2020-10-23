.. _label-linear-algebra-advanced:

Linear algebra
==============

We provide a set of useful function to perform basic linear algebra manipulations available by including ``castor/linalg.hpp``. The **linear algebra** library is mostly based on the well-known BLAS and LAPACK APIs. As such, it may be interfaced with any library respecting the same naming convention like ``openblas`` or ``mkl``. 

We want to bring the attention of the user on the fact that BLAS and LAPACK are originally Fortran libraries and consequently use a column ordering storage convention, unlike ``matrix`` which uses a row-major storage. While this behavior is fully transparent for the user as long as he is using the high-level **castor** API, he may find more details at :ref:`label-blaslapack-issue`.

We give below a few examples of use. More may be found in the ``demo/demo_linalg/`` subfolder. All linear algebra functions are described at :ref:`label-singular-eig-values-func`, :ref:`label-factorization-func` and :ref:`label-linear-solver-func`



Singular Value Decomposition
----------------------------
It is very easy to compute the SVD of a given matrix. The function ``svd`` simply returns a tuple containing the singular values ``S``, the left singular vectors ``U`` and finally the transposed right singular vectors ``Vt``.

.. code:: c++

    matrix<double> A({
        {1,1,-2,1,3,-1},
        {2,-1,1,2,1,-3},
        {1,3,-3,-1,2,1},
        {5,2,-1,-1,2,1},
        {-3,-1,2,3,1,3},
        {4,3,1,-6,-3,-2}
    });
    matrix<> S, U, Vt;
    std::tie(S,U,Vt) = svd(A,"vect");
    disp(norm(A - mtimes(U,mtimes(diag(S),Vt)),"inf"));

.. code:: text

    4.44089e-15

Note that it is equally easy to compute the (eventually generalized) eigenvalues and eigenvectors (see :ref:`label-eig`) or the rank-revealing SVD (see :ref:`label-qrsvd`).


Solving a linear system
---------------------------
We provide two functions for the resolution of a linear system of equations. ``linsolve`` performs the exact inversion using the LAPACK library while ``gmres`` computes the solution iteratively using the GMRES algorithm. Note that **both** have support for multiple right-hand-sides.

The use of ``linsolve`` is straightforward. The left- and right-hand-sides should be provided in full ``matrix`` form and the function returns the ``matrix`` containing the solutions.

.. code:: c++

    // left-hand-side
    matrix<double> A({
        {1,1,-2,1,3,-1},
        {2,-1,1,2,1,-3},
        {1,3,-3,-1,2,1},
        {5,2,-1,-1,2,1},
        {-3,-1,2,3,1,3},
        {4,3,1,-6,-3,-2}
    });
    // we change the original matrix for this demo
    for(auto i=0; i<6; ++i) A(i,i) = 20.; 
    // right-hand-side, must be a column
    matrix<double> B(6,1,{4,20,-15,-3,16,-27});
    // solving
    auto X = linsolve(A,B);
    disp(X);

.. code:: text

   -0.15882  
    0.83718  
   -0.93171  
   -0.29017  
    1.15147  
   -1.31156

Compared to  ``linsolve``, one needs to specifiy the tolerance and the maximum number of iterations. The user may also provide a preconditionner, use its own definition of the matrix-vector product, etc., as specified in the documentation of :ref:`label-gmres`.

.. code:: c++

    X = gmres(A,B,1e-3,6);
    disp(X);

.. code:: text

    Start GMRES using MGCR implementation (Multiple Generalized Conjugate Residual):
    + Iteration 1 in 3.5977e-05 seconds with relative residual 0.296391.
    + Iteration 2 in 2.8044e-05 seconds with relative residual 0.0545852.
    + Iteration 3 in 2.5446e-05 seconds with relative residual 0.0119836.
    + Iteration 4 in 3.3522e-05 seconds with relative residual 0.000506301.
    GMRES converged at iteration 4 to a solution with relative residual 0.000506301.
       -0.15823  
        0.83694  
       -0.93228  
       -0.29023  
        1.15188  
       -1.31119

**Remark:** Exact convergence is always achieved when the number of iterations reaches the dimension of the matrix.


Adaptive Cross Approximation
----------------------------
The ACA algorithm performs a low-rank approximation of rank ``k`` a given matrix ``C`` of size ``m x n`` under the form ``C = mtimes(A,tranpose(B))`` where ``A`` is ``m x k`` and ``B`` is ``n x k``. The rank ``k`` is automatically obtained following a user-defined accuracy parameter ``tol``. Note that if ``tol`` is chosen too small, using the ACA algorithm may be useless ...
**Warning** : our implementation directly returns ``Bt = tranpose(B)``.

.. code:: c++

    matrix<> A, Bt;
    matrix<> C = mtimes(rand(100,20),rand(20,50));
    std::tie(A,Bt) = aca(C,1e-3);
    disp(size(A),2);
    disp(size(Bt),2);
    disp(norm(C - mtimes(A,Bt),"inf"),2);

.. code:: text

    Matrix 1x2 of type 'm' (16 B):
    100   21  
    Matrix 1x2 of type 'm' (16 B):
    21  50  
    Object of type 'd':
    3.81917e-14

**Remark:** The result may vary depending on your random number generator.



.. _label-blaslapack-issue:

Overload of some BLAS and LAPACK functions
------------------------------------------

Some functions from the BLAS and LAPACK APIs have been directly interfaced using a similar naming convention. 

The BLAS-part is in fact a straightforward overlay over the C BLAS API which enables the possibility to use row-major ordering. For example, ``tgemm`` computes the matrix-matrix product using the C BLAS API (see :ref:`label-tgemm-blas`) and provides a unique interface to the four corresponding functions (each one corresponding to one of the four types: ``float``, ``std::complex<float>``, ``double``, ``std::complex<double>``). Its main purpose is to hide the use of raw pointers to the underlying data. For information, such raw pointers can be obtained as explained below.

.. code:: c++

    matrix<double> A({1,2,3,4});
    double *pA = &A(0); // returns the adress of the first element of A

However, the LAPACK-part is a *direct* overlay over the Fortran LAPACK API. Under the hood, the functions of **castor** convert the ``matrix<>`` to a ``std::vector<>`` where the data is stored with the right ordering thanks to the ``mat2lpk`` function (see :ref:`label-mat2lpk`). The result is then converted back to a ``matrix<>`` with the ``lpk2mat`` function (see :ref:`label-lpk2mat`). The main consequence is that the functions following a naming convention close to the LAPACK one (see for example :ref:`label-tgesdd`, :ref:`label-tgeqrf`, etc.) only accept ``std::vector<>`` as input.
