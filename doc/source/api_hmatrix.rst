
API
===

.. _label-aca-hmatrix:

aca
---
.. doxygenfunction:: aca(matrix<std::size_t> I, matrix<std::size_t> J, std::function<matrix<double>(matrix<std::size_t>, matrix<std::size_t>)> const &fct, double tol)
    :project: castor
.. doxygenfunction:: aca(matrix<std::size_t> I, matrix<std::size_t> J, std::function<matrix<std::complex<double>>(matrix<std::size_t>, matrix<std::size_t>)> const &fct, double tol)
    :project: castor


.. _label-full-hmatrix:

full
----
.. doxygenfunction:: full(hmatrix<T> const &Ah)
    :project: castor

.. _label-gmres-hmatrix:

gmres
-----
.. doxygenfunction:: gmres(hmatrix<T> const &Ah, matrix<T> const &B, double tol = 1e-6, std::size_t maxit = 10, hmatrix<T> const &Ahm1 = hmatrix<T>(), matrix<T> const &X0 = matrix<T>())
    :project: castor

.. _label-inv-hmatrix:

inv
---
.. doxygenfunction:: inv(hmatrix<T> const &Ah)
    :project: castor


.. _label-linsolve-hmatrix:

linsolve
--------
.. doxygenfunction:: linsolve(hmatrix<T> const &Ah, matrix<T> const &B)
    :project: castor


.. _label-lu-hmatrix:

lu
--
.. doxygenfunction:: lu(hmatrix<T> const &Ah)
    :project: castor


.. _label-mtimes-hmatrix:

mtimes
------
.. doxygenfunction:: mtimes(hmatrix<T> const &Ah, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: mtimes(matrix<T> const &A, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: mtimes(hmatrix<T> const &Ah, matrix<T> const &B)
    :project: castor


.. _label-rank2-hmatrix:

rank2
-----
.. doxygenfunction:: rank2(matrix<T> const &A, double tol)
    :project: castor


.. _label-size_hmatrix:

size
----
.. doxygenfunction:: size(hmatrix<T> const &Ah)
    :project: castor

See :ref:`label-length`, :ref:`label-numel`.


.. _label-spy-hmatrix:

spy
---
.. doxygenfunction:: spy(hmatrix<T> const &Ah)
    :project: castor


.. _label-tgeabm-hmatrix:

tgeabm
------
.. doxygenfunction:: tgeabm(T alpha, matrix<T> const &A, matrix<T> const &B, T beta, hmatrix<T> &Ch)
    :project: castor


.. _label-tgemm-hmatrix:

tgemm
-----
.. doxygenfunction:: tgemm(T alpha, hmatrix<T> const &Ah, hmatrix<T> const &Bh, T beta, hmatrix<T> &Ch)
    :project: castor


.. _label-transpose-hmatrix:

transpose
---------
.. doxygenfunction:: transpose(hmatrix<T> const &Ah)
    :project: castor


.. _label-operators-hmatrix:

Operators
+++++++++

.. _label-operator<<-hmatrix:

operator<<
----------
.. doxygenfunction:: operator<<(std::ostream &flux, hmatrix<T> const &Ah)
    :project: castor


.. _label-operator+-hmatrix:

operator+
---------
.. doxygenfunction:: operator+(T a, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator+(hmatrix<T> const &Ah, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator+(hmatrix<T> const &Ah, T b)
    :project: castor


.. _label-operator--hmatrix:

operator-
---------
.. doxygenfunction:: operator-(hmatrix<T> const &Ah)
    :project: castor
.. doxygenfunction:: operator-(T a, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator-(hmatrix<T> const &Ah, hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator-(hmatrix<T> const &Ah, T b)
    :project: castor


.. _label-operator*-hmatrix:

operator*
---------
.. doxygenfunction:: operator*(hmatrix<T> const &Ah, T b)
    :project: castor
.. doxygenfunction:: operator*(T a, hmatrix<T> const &Bh)
    :project: castor


.. _label-operator/-hmatrix:

operator/
---------
.. doxygenfunction:: operator/(hmatrix<T> const &Ah, T b)
    :project: castor


.. _label-operator+=-hmatrix:

operator+=
----------
.. doxygenfunction:: operator+=(hmatrix<T> const &Bh)
    :project: castor
.. doxygenfunction:: operator+=(T b)
    :project: castor



.. _label-operator*=-hmatrix:

operator*=
----------
.. doxygenfunction:: operator*=(T b)
    :project: castor
