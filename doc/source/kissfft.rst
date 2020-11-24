.. _label-using-kissfft:

Fast Fourier Transform
======================

The **castor** framework wraps the `kissfft library <https://github.com/mborgerding/kissfft>`_ in a standalone ``kissfft.hpp`` header file so it is not necessary to download it. Currently, only **single precision** arithmetic is supported. We describe below a basic example.

A minimum working file can be found below.

.. code:: c++

    #include "castor/matrix.hpp"
    #include "castor/kissfft.hpp"

    int main()
    {
        matrix<std::complex<float>> A = rand(3,4) + M_1I*rand(3,4);
        // your code here
        return 0;
    }

In the code above, we created a random complex single precision matrix with twelve elements. Since only the one-dimensional fft is supported, the matrix ``A`` will be treated as a one-dimensional array where the first row comes first, etc. Now we declare a second matrix ``B`` where we store the result of the forward-backward discrete Fourier transform.

.. code:: c++

    matrix<std::complex<float>> ifft(fft(A));

While ``A`` is a two-dimensional array, this is not the case for ``B``. In order to access the one-dimensional data of ``A``, we use a view on the 12 elements using ranged access (see :ref:`label-view`) and we compare it to ``B``.

.. code:: c++

    disp(norm(eval(A(range(0,12))) - B, "inf");

The output should be something more or less like

.. code:: text

    (9.31323e-08,0)

Please note that this code is part of the ``overview_kissfft.cpp`` file which can be found in the ``demo/demo_kissfft`` subfolder.