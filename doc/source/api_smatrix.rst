.. _label-smatrix-API:

API
===

.. _label-check-smatrix:

check
-----
.. doxygenfunction:: check(smatrix<T> &As)
    :project: castor

See :ref:`label-class-matrix`.


.. _label-clear-smatrix:

clear
-----
.. doxygenfunction:: clear(smatrix<T> &As)
    :project: castor

See :ref:`label-class-matrix`.


.. _label-disp-smatrix:

disp
----
.. doxygenfunction:: disp(smatrix<T> const &As, int info = 2, std::ostream &flux = std::cout, std::size_t r = 3)
    :project: castor

See :ref:`label-operator<<-smatrix`.


.. _label-find-smatrix:

find
----
.. doxygenfunction:: find(smatrix<T> const &As)
    :project: castor

See :ref:`label-index-smatrix`, :ref:`label-values-smatrix`.


.. _label-full-smatrix:

full
----
.. doxygenfunction:: full(smatrix<T> const &As)
    :project: castor

See :ref:`label-sparse-smatrix`.


.. _label-index-smatrix:

index
-----
.. doxygenfunction:: index(smatrix<T> const &As)
    :project: castor

.. _label-mtimes-smatrix:

mtimes
------
.. doxygenfunction:: mtimes(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: mtimes(smatrix<R> const &As, matrix<S> const &B)
    :project: castor

See :ref:`label-tgemm-naive`, :ref:`label-kron`.


.. _label-nnz-smatrix:

nnz
---
.. doxygenfunction:: nnz(smatrix<T> const &As)
    :project: castor

See :ref:`label-find`, :ref:`label-size-smatrix`.


.. _label-numel-smatrix:

numel
-----
.. doxygenfunction:: numel(smatrix<T> const &As)
    :project: castor

See :ref:`label-size-smatrix`, :ref:`label-nnz-smatrix`.


.. _label-reshape-smatrix:

reshape
-------
.. doxygenfunction:: reshape(smatrix<T> const &As, std::size_t m, std::size_t n)
    :project: castor

See :ref:`label-transpose-smatrix`.

.. _label-size-smatrix:

size
----
.. doxygenfunction:: size(smatrix<T> const &As, int dim)
    :project: castor
.. doxygenfunction:: size(smatrix<T> const &As)
    :project: castor

See :ref:`label-numel-smatrix`, :ref:`label-nnz-smatrix`.


.. _label-sparse-smatrix:

sparse
------

.. doxygenfunction:: sparse(matrix<std::size_t> const &L, matrix<T> const &V, std::size_t m, std::size_t n)
    :project: castor
.. doxygenfunction:: sparse(matrix<std::size_t> const &I, matrix<std::size_t> const &J, matrix<T> const &V)
    :project: castor
.. doxygenfunction:: sparse(matrix<std::size_t> const &I, matrix<std::size_t> const &J, matrix<T> const &V, std::size_t m, std::size_t n)
    :project: castor
.. doxygenfunction:: sparse(matrix<T> const &A)
    :project: castor

See :ref:`label-full-smatrix`.


.. _label-speye-smatrix:

speye
-----
.. doxygenfunction:: speye(matrix<std::size_t> const &S)
    :project: castor
.. doxygenfunction:: speye(std::size_t m, long n = -1)
    :project: castor

See :ref:`label-spzeros-smatrix`, :ref:`label-spones-smatrix`, :ref:`label-sprand-smatrix`, :ref:`label-eye`.


.. _label-spones-smatrix:

spones
------
.. doxygenfunction:: spones(matrix<std::size_t> const &S)
    :project: castor
.. doxygenfunction:: spones(std::size_t m, long n = -1)
    :project: castor

See :ref:`label-spzeros-smatrix`, :ref:`label-speye-smatrix`, :ref:`label-sprand-smatrix`, :ref:`label-ones`.


.. _label-sprand-smatrix:

sprand
------
.. doxygenfunction:: sprand(matrix<std::size_t> const &S, bool seed = false)
    :project: castor
.. doxygenfunction:: sprand(std::size_t m, long n = -1, bool seed = false)
    :project: castor

See :ref:`label-spzeros-smatrix`, :ref:`label-speye-smatrix`, :ref:`label-spones-smatrix`, :ref:`label-rand`.


.. _label-spzeros-smatrix:

spzeros
-------
.. doxygenfunction:: spzeros(std::size_t m, long n = -1)
    :project: castor

See :ref:`label-spones-smatrix`, :ref:`label-speye-smatrix`, :ref:`label-sprand-smatrix`, :ref:`label-zeros`.


.. _label-transpose-smatrix:

transpose
---------
.. doxygenfunction:: transpose(smatrix<T> const &As)
    :project: castor

See :ref:`label-reshape-smatrix`.


.. _label-values-smatrix:

values
------
.. doxygenfunction:: values(smatrix<T> const &As)
    :project: castor

See :ref:`label-index-smatrix`, :ref:`label-find-smatrix`.



.. _label-operator<<-smatrix:

operator<<
----------
.. doxygenfunction:: operator<<(std::ostream &flux, smatrix<T> const &As)
    :project: castor

See :ref:`label-disp-smatrix`.


.. _label-operator+-smatrix:

operator+
---------
.. doxygenfunction:: operator+(smatrix<R> const &As, S b)
    :project: castor
.. doxygenfunction:: operator+(R a, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator+(smatrix<R> const &As, matrix<S> const &B)
    :project: castor
.. doxygenfunction:: operator+(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator+(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor



.. _label-operator--smatrix:

operator-
---------
.. doxygenfunction:: operator-(smatrix<R> const &As, S b)
    :project: castor
.. doxygenfunction:: operator-(R a, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator-(smatrix<R> const &As, matrix<S> const &B)
    :project: castor
.. doxygenfunction:: operator-(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator-(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator-(smatrix<T> const &As)
    :project: castor



.. _label-operator*-smatrix:

operator*
---------
.. doxygenfunction:: operator*(smatrix<R> const &As, S Bs)
    :project: castor
.. doxygenfunction:: operator*(R As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator*(smatrix<R> const &As, matrix<S> const &B)
    :project: castor
.. doxygenfunction:: operator*(matrix<R> const &A, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator*(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor


.. _label-operator/-smatrix:

operator/
---------
.. doxygenfunction:: operator/(smatrix<R> const &As, S Bs)
    :project: castor
.. doxygenfunction:: operator/(matrix<R> const &As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator/(smatrix<R> const &As, matrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator/(R As, smatrix<S> const &Bs)
    :project: castor
.. doxygenfunction:: operator/(smatrix<R> const &As, smatrix<S> const &Bs)
    :project: castor

