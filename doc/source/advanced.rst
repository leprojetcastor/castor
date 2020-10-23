.. _label-advanced:

Advanced
========

We present here some advanced features of the ``matrix`` class.

Working with complex numbers
++++++++++++++++++++++++++++

It is possible to manipulate complex data with the ``matrix`` class. The code below creates a complex-double ``matrix`` using the predefined constant ``M_1I`` which is equal to ``std::complex<double>(0,1)`` and two matrices containing respectively the real and the imaginary part.


.. code:: c++

    matrix<> Re = {1,0,-1,0};
    matrix<> Im = {0,1,0,-1};

    disp("Imaginary number matrix ('1i') : ");
    disp(M_1I);

    disp("Complex matrix : ");
    auto Ac = Re + M_1I*Im;
    disp(Ac);


.. code:: text

  Imaginary number matrix ('1i') :
  (0,1)
  Complex matrix :
  (1.00000,0.00000)          (0.00000,1.00000)         (-1.00000,0.00000)         (0.00000,-1.00000)

We can compute the conjugate of the elements, the absolute value, or extract the imaginary part.

.. code:: c++

  disp(conj(Ac),2);
  disp(abs(Ac),2);
  disp(imag(Ac),2);

.. code:: text

  Matrix 1x4 of type 'St7complexIdE' (64 B):
        (1.00000,-0.00000)         (0.00000,-1.00000)        (-1.00000,-0.00000)          (0.00000,1.00000)  
  Matrix 1x4 of type 'd' (32 B):
      1.00000      1.00000      1.00000      1.00000  
  Matrix 1x4 of type 'd' (32 B):
            0      1.00000            0     -1.00000


Logical matrices
++++++++++++++++

Because of the special behaviour of ``std::vector<bool>``,  **matrix** does not use the native logical type ``bool`` of the C++ standard library, but uses instead ``std::uint8_t`` which corresponds to a ``unsigned short int``, or ``char``. Consequently, the ``bool`` values are converted to ``std::uint8_t`` where the value ``false`` is converted to the value ``0`` and ``true`` to the value ``1``. A logical ``matrix`` is displayed like:

.. code:: c++

    matrix<int> A = eye(3);
    matrix<int> B = ones(3);
    disp("A && B :");
    disp(A && B);

.. code:: text

    A && B :
    1  0  0
    0  1  0
    0  0  1

**WARNING (very important):** a ``matrix<logical>`` *remains intrinsically a* ``matrix<std::uint8_t>`` meaning that it behaves like one. This behavior is not natural and will be corrected in a future release.


Using some algorithms
+++++++++++++++++++++

The ``matrix`` class comes with a lot of built-in algorithms (see :ref:`label-algorithms` for a full list). First, we create some random integer data in the range ``[0,10[``.

.. code:: c++

  auto A = cast<int>(10*rand<>(1,10));
  disp(A,2,std::cout,10,10);

.. code:: text

  Matrix 1x10 of type 'i' (40 B):
  8  3  7  7  9  1  3  7  2  5

**Remark:** In the code above, we first call :ref:`label-rand` for the default ``double`` type, then we multiply by ``10`` and finally :ref:`label-cast` the result to integers. In fact, the :ref:`label-rand` function generates intrinsically numbers of type ``double`` in the range ``[0,1[`` which are converted afterward in the template type. Unfortunately, calling ``rand<int>`` will cast those numbers to zero. What we did is first generate ``double`` numbers in the range ``[0,10[``, then cast them to ``int``.

What are the minimum, average, standard deviation, and the maximum values of ``A`` ?

.. code:: c++

  disp(min(A));
  disp(mean<double>(A));
  disp(stddev<double>(A));
  disp(max(A));


.. code:: text

  1
  5.2
  2.63818
  9

**Remark:** Do not forget, in that particular case, to give the output type for the result of ``mean`` and ``stddev``. Otherwise, the result will be cast in the template type which is ``int``!

The :ref:`label-unique` elements of ``A`` are

.. code:: c++

  disp(unique(A),2,std::cout,10,10);

.. code:: text

  Matrix 1x7 of type 'i' (28 B):
  1  2  3  5  7  8  9

Now let us :ref:`label-sort` the elements.

.. code:: c++

  auto A_sorted = sort(A);
  disp(A_sorted,2,std::cout,10,10);

.. code:: text

  Matrix 1x10 of type 'i' (40 B):
  1  2  3  3  5  7  7  7  8  9  

We create a second smaller vector and we search the common elements.

.. code:: c++

  auto B = cast<int>(10*rand<>(1,4));
  disp(B,2);
  disp(intersect(A,B),2);

.. code:: text

  Matrix 1x4 of type 'i' (16 B):
  4  6  3  5  
  Matrix 1x2 of type 'i' (8 B):
  3  5  

We obtain the expected result.


View
++++

The **matrix** provides operators to extract a submatrix from a matrix instance or to assign values of a submatrix to a matrix. 
The operator ``()`` can take a list of indices to return an instance of ``class view``. 

As the ``class view`` contains a reference to the set of values corresponding to the list of indices, it is **necessary** to call the method ``eval`` to use the viewed submatrix.

It is possible to give the list of indices in linear indexing **A(L)** or bilinear indexing **A(I,J)**. The external functions all, col, range and row are useful to describe indices.

.. code:: c++

    matrix<> A = {{0, 1, 2, 4},
                  {4, 5, 6, 7},
                  {8, 9,10,11}};

    disp("Linear indexing");
    disp("  extracting :");
    matrix<> B = eval(A({1,3,5}));
    disp(B);

    disp("  assigning :"); 
    A(range(0,4)) = 0;
    disp(A);

    disp("Bilinear indexing");
    disp("  extracting :");
    B = eval(A(1,range(0,3)));;
    disp(B);

    disp("  assigning :"); 
    A({0,2}, col(A)) = 10;
    disp(A);

.. code:: text

  Linear indexing
    extracting :
      1.00000      4.00000      5.00000
    assigning :
            0            0            0            0
      4.00000      5.00000      6.00000      7.00000
      8.00000      9.00000     10.00000     11.00000
  Bilinear indexing
    extracting :
      4.00000      5.00000      6.00000
    assigning :
     10.00000     10.00000     10.00000     10.00000
      4.00000      5.00000      6.00000      7.00000
     10.00000     10.00000     10.00000     10.00000   


Internal matrix tools
+++++++++++++++++++++

In addition to constructors and operators, the matrix class has few member functions (see internal tools list in :ref:`class matrix description <label-class-matrix>`). The functions can be user-friendy for native C++ coding.

Notes:

- The method ``size``  differs from the external function ``size`` for the parameter dim=0. Indeed, in that case, the method ``size`` returns the total number of elements in the matrix while the external function ``size`` returns the two-element vector containing the number of rows and columns in the matrix.

.. code:: c++

    matrix<> A = {{1,2,3},{4,5,6}};
    disp("Method size :");
    disp(A.size(0));
    disp(A.size(1));
    disp(A.size(2));
    disp("External function size :");
    disp(size(A));
    disp(size(A,1));
    disp(size(A,2));

.. code:: text

   Method size :
   6
   2
   3
   External function size :
   2  3
   2
   3

- The method ``resize`` is more efficient than the external function ``resize`` since the methode makes the maniplulation in-place which avoids a copy.  
