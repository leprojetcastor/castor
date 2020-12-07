
.. _label-singular-eig-values-func:

Singular and eigen values
+++++++++++++++++++++++++

This section regroups all of the linear algebra functions concerning singular and eigen values which are currently implemented within **castor**. In order to access these functions, **the user must include** ``castor/linalg.hpp``. Two levels of interface are available. The high-level interface provides functions with simplified arguments while the low-level interface is much closer to the BLAS/LAPACK API.

The user will find here high-level linear algebra functions. Some examples of use may also be found at :ref:`label-linear-algebra-advanced`.

.. _label-eig:

eig
---
.. doxygenfunction:: eig(matrix<T> const &A, matrix<T> const &B)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<double> const &A, matrix<double> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<float> const &A, matrix<float> const &B, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<double>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<std::complex<float>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<double> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: eig(matrix<float> const &A, std::string typ)
   :project: castor

See :ref:`label-qr`, :ref:`label-svd`.

.. _label-qrsvd:

qrsvd
-----
.. doxygenfunction:: qrsvd(matrix<T> const &A, matrix<T> const &B, float tol)
   :project: castor

See :ref:`label-svd`, :ref:`label-qr`.

.. _label-rank:

rank
----

.. doxygenfunction:: rank()
   :project: castor

See :ref:`label-svd`, :ref:`label-aca`.


.. _label-svd:

svd
---
.. doxygenfunction:: svd(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: svd(matrix<std::complex<double>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: svd(matrix<double> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: svd(matrix<std::complex<float>> const &A, std::string typ)
   :project: castor
.. doxygenfunction:: svd(matrix<float> const &A, std::string typ)
   :project: castor

See :ref:`label-rank`, :ref:`label-eig`, :ref:`label-qr`, :ref:`label-qrsvd`, :ref:`label-aca`.
