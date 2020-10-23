
.. _label-sparse-matrix:

Sparse matrix
=============

This library implements a lightweight class ``smatrix`` for the manipulation of sparse matrices, using solely the ``castor::matrix`` class without other external libraries. Many :ref:`label-builders-smatrix`, :ref:`label-functions-smatrix` and :ref:`label-operators-smatrix` are available in a similar fashion as for the :ref:`label-class-matrix` class. Most of the manipulations should be transparent to the user and the feeling should be really `Matlab-like <https://www.mathworks.com>`_. Basic examples are given in the corresponding :ref:`label-examples-smatrix` section.

The ``smatrix`` class is installed automatically with the other headers of **castor** library.

**This is work-in-progress** 

.. _label-examples-smatrix:

Examples
........

In order to use a ``smatrix``, the corresponding header ``castor/smatrix.hpp`` needs to be included. A minimum working example would thus look like this:

.. code:: c++

    #include "castor/matrix.hpp"
    #include "castor/smatrix.hpp"

    int main()
    {
        smatrix<> s;
        // your code here
        return 0;
    }

Building a smatrix
++++++++++++++++++

A ``smatrix`` will mostly behave like a normal full ``matrix``. For example, let us create a ``4 x 4`` ``smatrix`` of ``double``.

.. code:: c++

    smatrix<> As;
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 0 elements (0 B):
    -empty-

The newly created ``smatrix`` is empty. It is rather easy to set the values of the entries. For example:

.. code:: c++

    As(1,2) = -0.5;
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 1 elements (8 B):
    (1,2)  -0.5

Now ``As`` contains one non-zero element with the bilinear index ``(1,2)`` (second line, third column). However, this way of filling a ``smatrix`` **is definitely not recommended** as it involves a lot of memory management and may affect dramatically the performances. One way to do so is to first create the matrix in *triplet* format ``(I,J,VALUE)`` and only afterward build the corresponding ``smatrix`` as illustrated below:

.. code:: c++

    matrix<std::size_t> I = {0,2,3};
    matrix<std::size_t> J = {1,1,2};
    matrix<> V = {0.5, -1/3., M_PI};

    As = smatrix<>(I,J,V,4,4);
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 3 elements (24 B):
    (0,1)  0.5
    (2,1)  -0.333333
    (3,2)  3.14159

Yes, we reaffected ``As`` to a new ``smatrix``. The old data is automatically discarded so one should be careful when performing such an action. As for ``matrix``, it is possible to :ref:`label-clear-smatrix` the content of a ``smatrix`` (the object is reinitialized). Let us now add a ``0``.

.. code:: c++

    As(3,3) = 0.; // but, why ?
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 4 elements (32 B):
    (0,1)  0.5
    (2,1)  -0.333333
    (3,2)  3.14159
    (3,3)  0

A zero value is added. In order to clean a ``smatrix``, a simple call to ``check`` is sufficient:

.. code:: c++

    check(As);
    disp(As,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 3 elements (24 B):
    (0,1)  0.5
    (2,1)  -0.333333
    (3,2)  3.14159

Everything went back to normal! Now, let us use one of the :ref:`label-builders-smatrix` in order to create an identity sparse matrix. It is also possible to convert back to the *triplet* format.

.. code:: c++

    auto Bs = speye<>(4,4);
    disp(Bs,2);
    matrix<std::size_t> IB,JB;
    matrix<> VB;
    std::tie(IB,JB,VB) = find(Bs);
    disp(transpose(vertcat(vertcat(IB,JB),VB)),2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 4 elements (32 B):
    (0,0)  1
    (1,1)  1
    (2,2)  1
    (3,3)  1
    Matrix 4x3 of type 'd' (96 B):
              0            0      1.00000  
        1.00000      1.00000      1.00000  
        2.00000      2.00000      1.00000  
        3.00000      3.00000      1.00000 

The matrices ``IB,JB,VB`` are returned as *line* vectors. To obtain a better display, we concatenated them vertically and tranposed the result.


Basic manipulations
+++++++++++++++++++

In this section, we start with start from scratch so everything written in the previous section should be discarded from your ``main`` function. Let us create two matrices

.. code:: c++

    smatrix<> As = speye(4,4);
    
    matrix<std::size_t> I({1,1,2,2}), J({1,2,1,2});
    matrix<> V({1.,1.,1.,1.});
    smatrix<> Bs = smatrix<>(I,J,V,4,4);

``As`` is a ``4 x 4`` identity matrix and ``Bs`` is a matrix with the interior filled with ones. Here is an example of basic manipulations:

.. code:: c++

    auto Cs = 1.5*As - Bs/2.;
    disp(Cs,2);

.. code:: text

    Sparse matrix 4x4 of type 'd' with 6 elements (48 B):
    (0,0)  1.5
    (1,1)  1
    (1,2)  -0.5
    (2,1)  -0.5
    (2,2)  1
    (3,3)  1.5

What is the number of non-zero elements, again ?

.. code:: c++

    std::cout << "nnz(Cs) = " << nnz(Cs) << std::endl;

.. code:: text

    nnz(Cs) = 6

It is possible to get the value of any entry:

.. code:: c++ 

    std::cout << "Cs(1,2) = " << Cs(1,2) << std::endl;
    std::cout << "Cs(1,3) = " << Cs(1,3) << std::endl;

.. code:: text

    Cs(1,2) = -0.5
    Cs(1,3) = 0

Now, let us multiply ``Cs`` by a ``4 x 4`` full ``matrix``:

.. code:: c++

    auto D = mtimes(Cs,ones<>(4));
    disp(D,2);  // :)

.. code:: text

    Matrix 4x4 of type 'd' (128 B):
        1.50000      1.50000      1.50000      1.50000  
        0.50000      0.50000      0.50000      0.50000  
        0.50000      0.50000      0.50000      0.50000  
        1.50000      1.50000      1.50000      1.50000

One last manipulation and we are good for this example.

.. code:: c++

    Cs(0,3) = M_PI;
    auto Es = Cs - transpose(Cs);
    check(Es); // drop the zeros
    disp(Es,2);
    disp(full(Es,2));

.. code:: text

    Sparse matrix 4x4 of type 'd' with 2 elements (16 B):
    (0,3)  3.14159
    (3,0)  -3.14159
    Matrix 4x4 of type 'd' (128 B):
              0            0            0      3.14159  
              0            0            0            0  
              0            0            0            0  
       -3.14159            0            0            0
