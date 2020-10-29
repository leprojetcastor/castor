Getting started
===============

This section gives an example of use of **matrix** library and the procedure to compile it.

Example
-------

.. code::

    #include <iostream>
    #include "castor/matrix.hpp"

    using namespace castor;

    int main (int argc, char* argv[])
    {
        matrix<float> A = {{ 1.0,  2.0,  3.0,  4.0},
                           { 5.0,  6.0,  7.0,  8.0},
                           { 9.0, 10.0, 11.0, 12.0}};
    
        matrix<double> B = eye(3,4);

        auto C = A + B;

        disp(C);
    
        return 0;
    }

This example displays the sum of two matrices with implicit cast :

.. code:: text

    2.00000      2.00000      3.00000      4.00000
    5.00000      7.00000      7.00000      8.00000
    9.00000     10.00000     12.00000     12.00000


.. _label-compilation:

Compilation 
-----------

Command line
++++++++++++

The library **matrix** is a header-only library. That actually means to compile program using **matrix**, you just have to specify to the compiler the path to the directory containing the headers. For instance, with GCC, the command to compile the above example (assuming it is contained in a file named ``main.cpp``) : 

.. code:: bash

    g++ -std=c++14 -I /path/to/castor/folder main.cpp -o main

IDE
+++

Simply enter the path to the directory containing the ``castor`` folder in your favorite IDE (refer to your IDE documentation).
