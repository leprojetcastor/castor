.. _label-basic:

Basics
======

In this section, the user will find informations about basic manipulations on the ``matrix`` object : creation (and destruction), etc. ``matrix`` is defined within the namespace ``castor`` meaning that in order to call the functions, their name should be preceeded by ``castor::``. In order to use directly their name, the user should add 

.. code:: c++

  using namespace castor;

in the preamble of the ``.cpp`` file. A minimum ``main`` file would then look like this:

.. code:: c++

  #include "castor/matrix.hpp"

  using namespace castor;

  int main()
  {
    matrix<double> A;
    // Write your code here
    return 0;
  }

For some more advanced features, see the corresponding section :ref:`label-advanced`. The user may also refer to the examples in the ``demo/*`` subdirectories of the main directory of the **castor** project.


Matrix creation (and destruction)
---------------------------------

There are two main paths to create a ``matrix`` object. The first way is to call one of the constructors explicitly. If no value is specified, the object will be filled with zeros. The second way it to call a *builder* (see :ref:`label-basic-builder` and :ref:`label-basic-display`).

From matrix constructor
+++++++++++++++++++++++

**Remark:** In the following, we will use a lot the :ref:`label-disp` function which is meant to produce a formatted output of the content of a ``matrix`` object.

The code below creates an empty ``matrix`` of type ``int``.

.. code:: c++

  matrix<int> A;
  disp(A,2);

.. code:: text

  Matrix 0x0 of type 'i' (0 B):
  -empty-

**Remark:** If the ``using namespace castor`` clause was not added in the preamble as explained at the beginning of this section, the code above becomes

.. code:: c++

  castor::matrix<int> A;
  castor::disp(A,2);
  // etc.

By passing the value "2" as a second argument to ``disp``, we can see that it is a ``matrix`` of size ``0x0`` of integer type whose data size is 0 Byte.

By passing as argument a single value of type *T*, a singleton ``matrix`` is created.

.. code:: c++

  matrix<int> A(M_PI);
  disp(A,2);

.. code:: text

  Matrix 1x1 of type 'i' (4 B):
  3

Here, ``A`` has been declared as an ``matrix`` of integers but a ``double`` containing the value of pi was passed as argument. As a consequence, it was cast to an ``int``, thus the value 3. Note that 4 Bytes is the size of an integer in C++ when standard compilation options are used.

Next, we initialize a matrix using an initialization-list. By passing a single list of values as argument to the constructor, a line-``matrix`` is created. By passing a list of a list, a ``matrix`` is created whose number of lines is the number of elements in the outer list and the number of columns is the number of elements in the inner lists. Please be careful that the number of elements in the inner list should be the same for all. These two options are illustrated below.

.. code:: c++

  matrix<float> A({1,2,3,4,5});       // matrix of floats
  matrix<>      B({{1,2,3},{4,5,6}}); // matrix of doubles
  disp(A,2);
  disp(B,2);

.. code:: text

  Matrix 1x5 of type 'f' (20 B):
      1.00000      2.00000      3.00000      4.00000      5.00000  
  Matrix 2x3 of type 'd' (48 B):
      1.00000      2.00000      3.00000  
      4.00000      5.00000      6.00000


Finally, it is possible to create a ``matrix`` by giving its dimensions and a fill-value. By default, the matrix is filled with 0s. In the example below, we create a ``2x3`` matrix filled with the value 4, then we modify one of the entries.

.. code:: c++

  matrix<> A(2,3,4.);
  disp(A,2);
  A(1,2) = -0.5;
  disp(A,2);

.. code:: text

  Matrix 2x3 of type 'd' (48 B):
      4.00000      4.00000      4.00000  
      4.00000      4.00000      4.00000  
  Matrix 2x3 of type 'd' (48 B):
      4.00000      4.00000      4.00000  
      4.00000      4.00000     -0.50000


Please refer to the constructors list in the :ref:`class matrix description <label-class-matrix>`. 


.. _label-basic-builder:

From builder
++++++++++++

We describe now some of the useful builders for the ``matrix`` class.

The code below creates a ``2x3`` matrix of ``double`` filled with zeros.

.. code:: c++

    matrix<> A = zeros(2,3);
    disp(A);

.. code:: text

   1.0000  1.0000  1.0000  
   1.0000  1.0000  1.0000  

**Remark:** this is equivalent to calling explicitly the ``matrix`` constructor.

The code below creates a ``1x10`` matrix of ``double`` initialized with linear spaced values :

.. code:: c++

    matrix<> A = linspace(0,1,10);
    disp(A,2);

.. code:: text

    Matrix 1x10 of type 'd' (80 B):
              0      0.11111      0.22222  ...      0.77778      0.88889      1.00000

For a ``2x2`` random matrix of ``float``, use

.. code:: c++

    matrix<float> A = rand<float>(2);
    disp(A,2);

.. code:: text

    Matrix 2x2 of type 'f' (16 B):
        0.84019      0.39438  
        0.78310      0.79844

This last result may differ depending on your random number generator.

Notes : 

- Matrices and vectors are objects of the matrix template class. A vector is considered as a (1xn) size by default or (nx1). 
- The template argument of class matrix is double by default. It is possible to specify type both for matrix constructors and builders.


Clear a matrix
++++++++++++++

If for some reason the content of a ``matrix`` needs to be cleared (for example, free some RAM), there are to possibilities. The first solution (the *clean one*) is to call the :ref:`label-clear` function.

.. code:: c++

  auto A=rand(5);
  disp(A,2);
  clear(A);
  disp(A,2);

.. code:: text

  Matrix 3x3 of type 'd' (72 B):
      0.84019      0.39438      0.78310  
      0.79844      0.91165      0.19755  
      0.33522      0.76823      0.27777  
  Matrix 0x0 of type 'd' (0 B):
  -empty-

The second one is to assign an empty ``matrix`` in place of an existing one.

.. code:: c++

  auto A=rand(3);
  disp(A,2);
  A = {};
  disp(A,2);

.. code:: text

  Matrix 3x3 of type 'd' (72 B):
      0.84019      0.39438      0.78310  
      0.79844      0.91165      0.19755  
      0.33522      0.76823      0.27777  
  Matrix 0x0 of type 'd' (0 B):
  -empty-



.. _label-basic-display:

Display
-------

A very useful function is the :ref:`label-disp` function. It produces a formatted output of the content of a ``matrix`` object with additional informations. Let us create a ``2x2`` random ``matrix``.

.. code:: c++

  auto A = rand<>(2);

The simplest call to :ref:`label-disp` displays the raw content without additional informations.

.. code:: c++

  disp(A);

.. code:: text

      0.84019      0.39438  
      0.78310      0.79844

In fact, this is equivalent to calling ``disp(A,1)``. The second (optional) argument determines the level of informations to be displayed. ``disp(A,0)`` will produce the same output as ``disp(A)`` but no end-of-line character is added to the output. ``disp(A,2)`` will add informations on the size, the type and the RAM storage of the ``matrix``, as illustrated before.

The third argument to :ref:`label-disp` is the output stream (``std::ostream``) which by default is the standard output ``std::cout``. Finally, the user may specify two additional arguments which are the number of lines and columns which need to be displayed. By default, their value is 3 meaning that the first 3 and last 3 element of each direction are displayed.

.. code:: c++

  auto A = rand(10);
  disp(A,2);

.. code:: text

  Matrix 10x10 of type 'd' (800 B):
      0.84019      0.39438      0.78310  ...      0.76823      0.27777      0.55397  
      0.47740      0.62887      0.36478  ...      0.71730      0.14160      0.60697  
      0.01630      0.24289      0.13723  ...      0.10881      0.99892      0.21826  
  ...
      0.53161      0.03928      0.43764  ...      0.73853      0.63998      0.35405  
      0.68786      0.16597      0.44010  ...      0.89337      0.35036      0.68667  
      0.95647      0.58864      0.65730  ...      0.81477      0.68422      0.91097

Now we modifiy a little bit the format.

.. code:: c++
  
  disp(A,2,std::cout,4,2);

.. code:: text

  Matrix 10x10 of type 'd' (800 B):
      0.84019      0.39438  ...      0.27777      0.55397  
      0.47740      0.62887  ...      0.14160      0.60697  
      0.01630      0.24289  ...      0.99892      0.21826  
      0.51293      0.83911  ...      0.29252      0.77136  
  ...
      0.23828      0.97063  ...      0.51254      0.66772  
      0.53161      0.03928  ...      0.63998      0.35405  
      0.68786      0.16597  ...      0.35036      0.68667  
      0.95647      0.58864  ...      0.68422      0.91097

**Note :** :ref:`label-disp` can also display the content of other variables :

.. code:: c++

  disp("pi value is :");
  disp(M_PI);

.. code:: text

    pi value is
    3.14159


Size and indexing 
-----------------

We describe now a few useful functions to begin manipulating matrices.

Size, length, numel
+++++++++++++++++++

The **size** function returns the two-element vector containing the number of rows and columns in the matrix. The result is *also* a ``matrix``.

.. code:: c++

  matrix<> A = eye(3,4);
  disp(size(A),2)

.. code:: text

  Matrix 1x2 of type 'm' (16 B):
  3  4

If a dimension is specified, :ref:`label-size` returns only the length of the specified dimensions :

.. code:: c++

  matrix<> A = eye(3,4);
  disp(size(A, 1));
  disp(size(A, 2));

.. code:: text

   3
   4

The :ref:`label-length` and :ref:`label-numel` functions returns respectively the maximum length and the number of elements in the matrix :

.. code:: c++

  matrix<> A = eye(3,4);
  disp(length(A));
  disp(numel(A));

.. code:: text

   4
   12


Accessing elements
++++++++++++++++++

The elements of a ``matrix`` can be accessed using either *linear* or *bilinear* indexing. 

*Linear* indexing consists in accessing the n-th element of the ``matrix`` in the order the data is stored. Since ``matrix`` uses a row-major layout, the rows of the ``matrix`` are concatenated one after the other. 

.. code:: c++

  matrix<> A = {{1,2,3},
                {4,5,6},
                {7,8,9}};
  disp(A(5));

.. code:: text

  6

The 5-th element of ``A`` thus holds the value 6 (the 6-th holds 7, etc.). Linear indexing is particularly useful when accessing the elements of a one-dimensional ``matrix`` (a *vector*).

**Remark:** The index numbering follows the C/C++ convention meaning that the indexes start a ``0`` and ends at ``n-1`` where ``n`` would be a dimension of the ``matrix`` (see :ref:`label-size`).

*Bilinear* indexing is the natural way to access the elements of a ``matrix``.

.. code:: c++

  disp(A(1,2));

.. code:: text

  6


Help
----

The function :ref:`label-help` allows you to display the documentation of a function at runtime. You have to give the complete path to the header file ``matrix.hpp`` in the ``documentationFiles`` variable. 

.. code:: c++

  documentationFiles =
  {
      "/complete/path/to/matrix.hpp"
  };

  help("size");

.. code:: text

  ============================ DOCUMENTATION ============================
  Help on "size":
  Size of array.

  S = size(A) for m-by-n matrix A returns the two-element vector [m,n]
  containing the number of rows and columns in the matrix.

  S = size(A,dim) returns the lengths of the specified dimensions dim.

  Example(s):
     matrix<> A = {{1,2,3},{4,5,6}};
     disp(size(A));
     disp(size(A,1));
     disp(size(A,2));

  See also:
    length, numel.
  =======================================================================


Basic operations
----------------

The ``matrix`` class is designed to be as easy of use as Matlab or Numpy arrays. As a consequence, many operators have been overloaded. We will describe here some basic manipulations with few of all of the available operators. They can be discovered at :ref:`label-operators`.

First, we create two matrices ``A`` and ``B``, we multiply the first one by ``M_PI`` and we add them.

.. code:: c++

  // create two 'double' matrices
  auto A = eye(2);
  auto B = eye(2);
  disp(A,2);
  disp(B,2);

  A *= M_PI;
  disp(A,2);

  auto C = A + B;
  disp(C,2);

.. code:: text

  Matrix 2x2 of type 'd' (32 B):
      1.00000            0  
            0      1.00000  
  Matrix 2x2 of type 'd' (32 B):
      1.00000            0  
            0      1.00000  
  Matrix 2x2 of type 'd' (32 B):
      3.14159            0  
            0      3.14159  
  Matrix 2x2 of type 'd' (32 B):
      4.14159            0  
            0      4.14159
      
Then, we create an orthogonal matrix ``L`` and we compute ``D = C - L*L'*C``. Orthogonal matrices are such that the matrix product with their transpose yields the identity matrix. Therefore, the result should be a null-matrix.

.. code:: c++

  double theta = 0.2;
  matrix<> L = {
    {std::cos(theta),-std::sin(theta)},
    {std::sin(theta),std::cos(theta)}
  };

  auto D = C - mtimes(L,mtimes(transpose(L),C));
  disp(D,2);

.. code:: text

  Matrix 2x2 of type 'd' (32 B):
          0            0  
          0            0  


We obtain the expected result. **Note that the matrix-matrix product is computed using** :ref:`label-mtimes`. Indeed, the ``*`` operator *does not* compute the matrix-matrix product but a term-by-term product. 

.. code:: c++

  auto A = ones(2);
  auto B = matrix<>(2,2,2.);
  disp(A*B,2);

.. code:: text

  Matrix 2x2 of type 'd' (32 B):
      2.00000      2.00000  
      2.00000      2.00000
