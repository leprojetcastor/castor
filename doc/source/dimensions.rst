.. _label-matrix-dimensions:

Matrix dimensions
+++++++++++++++++

These functions allow to recover the dimensions of a matrix : total number of elements (see :ref:`label-numel`), dimensions (see :ref:`label-size`), etc.

.. _label-length:

length
------
.. doxygenfunction:: length(matrix<T> const &A)
   :project: castor

See :ref:`label-numel`, :ref:`label-size`.

.. _label-nnz:

nnz
---
.. doxygenfunction:: nnz(matrix<T> const &A)
   :project: castor

See :ref:`label-find`, :ref:`label-size`.

.. _label-numel:

numel
-----
.. doxygenfunction:: numel(matrix<T> const &A)
   :project: castor

See :ref:`label-length`, :ref:`label-size`.

.. _label-size:

size
----
.. doxygenfunction:: size(matrix<T> const &A, int dim)
   :project: castor

.. doxygenfunction:: size(matrix<T> const &A)
   :project: castor

See :ref:`label-length`, :ref:`label-numel`.
