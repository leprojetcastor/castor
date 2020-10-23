
.. _label-factorization-func:

Factorization
+++++++++++++

This section regroups all of the linear algebra functions concerning factorization which are currently implemented within **castor**. In order to access these functions, **the user must include** ``castor/linalg.hpp``. Two levels of interface are available. The high-level interface provides functions with simplified arguments while the low-level interface is much closer to the BLAS/LAPACK API.

Some examples of use may also be found at :ref:`label-linear-algebra-advanced`.


.. _label-aca:

aca
---
.. doxygenfunction:: aca(matrix<std::complex<double>> const &M, double tol)
   :project: castor
.. doxygenfunction:: aca(matrix<double> const &M, double tol)
   :project: castor

See :ref:`label-svd`, :ref:`label-rank`.


.. _label-lu:

lu
--
.. doxygenfunction:: lu(matrix<std::complex<double>> const &A)
   :project: castor
.. doxygenfunction:: lu(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: lu(matrix<std::complex<float>> const &A)
   :project: castor
.. doxygenfunction:: lu(matrix<float> const &A)
   :project: castor

See :ref:`label-qr`, :ref:`label-linsolve`, :ref:`label-inv`.

.. _label-qr:

qr
--
.. doxygenfunction:: qr(matrix<std::complex<double>> const &A)
   :project: castor
.. doxygenfunction:: qr(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: qr(matrix<std::complex<float>> const &A)
   :project: castor
.. doxygenfunction:: qr(matrix<float> const &A)
   :project: castor

See :ref:`label-eig`, :ref:`label-svd`, :ref:`label-lu`.
