.. _label-matrix-manipulations:

Matrix manipulations
++++++++++++++++++++

Many functions are provided for the manipulation of matrices. For instance, it is possible to :ref:`label-cast` a ``matrix`` object into a different type, to concatenate matrices in all dimensions (see :ref:`label-cat`, etc.), :ref:`label-find` the non-zero elements, or to :ref:`label-transpose` it.

.. _label-all:

all
---
.. doxygenfunction:: all(matrix<T> const &A)
   :project: castor

See :ref:`label-row`, :ref:`label-col`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`, :ref:`label-cview`.


.. _label-cast:

cast
----
.. doxygenfunction:: cast
   :project: castor

See :ref:`label-class-matrix`.


.. _label-cat:

cat
---
.. doxygenfunction:: cat(int dim, R const &A, S const &B)
   :project: castor
.. doxygenfunction:: cat(int dim, R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: cat(int dim, matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: cat(int dim, matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-vertcat`, :ref:`label-horzcat`.


.. _label-clear:

clear
-----
.. doxygenfunction:: clear(matrix<T> &A)
   :project: castor

See :ref:`label-class-matrix`, :ref:`label-resize`.


.. _label-col:

col
---
.. doxygenfunction:: col(matrix<T> const &A)
   :project: castor

See :ref:`label-row`, :ref:`label-all`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`, :ref:`label-cview`.


.. _label-find:

find
----
.. doxygenfunction:: find(matrix<T> const &A)
   :project: castor

See :ref:`label-nnz`, :ref:`label-ind2sub`.


.. _label-get:

get
---
.. doxygenfunction:: get(matrix<T> const &A, matrix<std::size_t> const &I, matrix<std::size_t> const &J)
   :project: castor
.. doxygenfunction:: get(matrix<T> const &A, matrix<std::size_t> const &L)
   :project: castor

See :ref:`label-set`, :ref:`label-all`, :ref:`label-row`, :ref:`label-col`, :ref:`label-view`, :ref:`label-cview`.


.. _label-horzcat:

horzcat
-------
.. doxygenfunction:: horzcat
   :project: castor

See :ref:`label-vertcat`, :ref:`label-cat`.


.. _label-ind2sub:

ind2sub
-------
.. doxygenfunction:: ind2sub
   :project: castor

See :ref:`label-sub2ind`, :ref:`label-find`.


.. _label-isempty:

isempty
-------
.. doxygenfunction:: isempty
   :project: castor

See :ref:`label-isequal`, :ref:`label-isvector`.


.. _label-isequal:

isequal
-------
.. doxygenfunction:: isequal
   :project: castor

See :ref:`label-isempty`, :ref:`label-isvector`.


.. _label-isvector:

isvector
--------
.. doxygenfunction:: isvector
   :project: castor

See :ref:`label-isequal`, :ref:`label-isempty`.


.. _label-resize:

resize
------
.. doxygenfunction:: resize(std::size_t m, std::size_t n, T v = (T)NAN)
   :project: castor

See :ref:`label-reshape`.


.. _label-reshape:

reshape
-------
.. doxygenfunction:: reshape(std::size_t m, std::size_t n)
   :project: castor

See :ref:`label-resize`, :ref:`label-transpose`.


.. _label-row:

row
---
.. doxygenfunction:: row(matrix<T> const &A)
   :project: castor

See :ref:`label-all`, :ref:`label-col`, :ref:`label-get`, :ref:`label-set`, :ref:`label-view`, :ref:`label-cview`.


.. _label-set:

set
---
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &I, matrix<std::size_t> const &J, U b)
   :project: castor
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &L, matrix<U> const &B)
   :project: castor
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &I, matrix<std::size_t> const &J, matrix<U> const &B)
   :project: castor
.. doxygenfunction:: set(matrix<T> &A, matrix<std::size_t> const &L, U b)
   :project: castor

See :ref:`label-get`, :ref:`label-all`, :ref:`label-row`, :ref:`label-col`, ::ref:`label-view`, :ref:`label-cview`.


.. _label-sub2ind:

sub2ind
-------
.. doxygenfunction:: sub2ind
   :project: castor

See :ref:`label-ind2sub`, :ref:`label-find`.


.. _label-transpose:

transpose
---------
.. doxygenfunction:: transpose(matrix<T> const &A)
   :project: castor

See :ref:`label-reshape`.


.. _label-vertcat:

vertcat
-------
.. doxygenfunction:: vertcat
   :project: castor

See :ref:`label-horzcat`, :ref:`label-cat`.
