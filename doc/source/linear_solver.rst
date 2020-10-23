
.. _label-linear-solver-func:

Linear solver 
+++++++++++++

This section regroups all of the linear algebra functions concerning linear solver which are currently implemented within **castor**. In order to access these functions, **the user must include** ``castor/linalg.hpp``. Two levels of interface are available. The high-level interface provides functions with simplified arguments while the low-level interface is much closer to the BLAS/LAPACK API.

The user will find here high-level linear algebra functions. Some examples of use may also be found at :ref:`label-linear-algebra-advanced`.


.. _label-inv:

inv
---
.. doxygenfunction:: inv(matrix<T> const &A)
   :project: castor

See :ref:`label-linsolve`, :ref:`label-pinv`, :ref:`label-gmres`.


.. _label-linsolve:

linsolve
--------
.. doxygenfunction:: linsolve(matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B)
   :project: castor
.. doxygenfunction:: linsolve(matrix<double> const &A, matrix<double> const &B)
   :project: castor
.. doxygenfunction:: linsolve(matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B)
   :project: castor
.. doxygenfunction:: linsolve(matrix<float> const &A, matrix<float> const &B)
   :project: castor

See :ref:`label-inv`, :ref:`label-pinv`, :ref:`label-gmres`.


.. _label-pinv:

pinv
----
.. doxygenfunction:: pinv(matrix<T> const &A)
   :project: castor

See :ref:`label-inv`, :ref:`label-pinv`.
