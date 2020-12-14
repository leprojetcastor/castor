.. _label-algorithms:

Algorithms
++++++++++

Standard algorithms over data stored in matrices may be found here. Among *many* others, it is possible to :ref:`label-sort` a ``matrix``, find the :ref:`label-unique` elements, the :ref:`label-median` of the data, perform the :ref:`label-dot` product between two matrices, find the maximum value with :ref:`label-max`, etc.


.. _label-argintersect:

argintersect
------------
.. doxygenfunction:: argintersect
   :project: castor

See :ref:`label-intersect`, :ref:`label-argsetdiff`, :ref:`label-argunique`.

.. _label-argmax:

argmax
------
.. doxygenfunction:: argmax(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: argmax(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-maximum`, :ref:`label-argmin`, :ref:`label-argsort`.

.. _label-argmin:

argmin
------
.. doxygenfunction:: argmin(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: argmin(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-min`, :ref:`label-minimum`, :ref:`label-argmax`, :ref:`label-argsort`.

.. _label-argsetdiff:

argsetdiff
----------
.. doxygenfunction:: argsetdiff
   :project: castor

See :ref:`label-setdiff`, :ref:`label-argintersect`, :ref:`label-argunique`.

.. _label-argsort:

argsort
-------
.. doxygenfunction:: argsort
   :project: castor

See :ref:`label-sort`, :ref:`label-argmin`, :ref:`label-argmax`.

.. _label-argunique:

argunique
---------
.. doxygenfunction:: argunique
   :project: castor

See :ref:`label-unique`, :ref:`label-argsetdiff`, :ref:`label-argintersect`.

.. _label-conv:

conv
----
.. doxygenfunction:: conv
   :project: castor

See :ref:`label-fftconv`.

.. _label-cross:

cross
-----
.. doxygenfunction:: cross
   :project: castor

See :ref:`label-dot`.

.. _label-cumprod:

cumprod
-------
.. doxygenfunction:: cumprod(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: cumprod(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-cumsum`, :ref:`label-prod`.

.. _label-cumsum:

cumsum
------
.. doxygenfunction:: cumsum(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: cumsum(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-cumprod`, :ref:`label-sum`.

.. _label-diff:

diff
----
.. doxygenfunction:: diff(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: diff(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-sum`, :ref:`label-prod`.

.. _label-dot:

dot
---
.. doxygenfunction:: dot
   :project: castor

See :ref:`label-cross`.

.. _label-fftconv:

fftconv
-------
.. doxygenfunction:: fftconv
   :project: castor

See :ref:`label-conv`, :ref:`label-fft`.


.. _label-gmres:

gmres
-----
.. doxygenfunction::  gmres(matrix<T> const &A, matrix<T> const &B, double tol = 1e-6, std::size_t maxit = 10, std::function<matrix<T>(matrix<T> const&)> const &Am1 = std::function<matrix<T>(matrix<T> const&)>(), matrix<T> const &X0 = matrix<T>())
   :project: castor
.. doxygenfunction::  gmres(matrix<T> const &A, matrix<T> const &B, double tol, std::size_t maxit, matrix<T> const &Am1, matrix<T> const &X0 = matrix<T>())
   :project: castor
.. doxygenfunction:: gmres(std::function<matrix<T>(matrix<T> const&)> const &A, matrix<T> const &B, double tol = 1e-6, std::size_t maxit = 10, std::function<matrix<T>(matrix<T> const&)> const &Am1 = std::function<matrix<T>(matrix<T> const&)>(), matrix<T> const &X0 = matrix<T>())
   :project: castor



See :ref:`label-linsolve`.

.. _label-intersect:

intersect
---------
.. doxygenfunction:: intersect
   :project: castor

See :ref:`label-argintersect`, :ref:`label-setdiff`, :ref:`label-unique`, :ref:`label-union2`.

.. _label-kron:

kron
----
.. doxygenfunction:: kron(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: kron(matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: kron(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-mtimes`.

.. _label-max:

max
---
.. doxygenfunction:: max(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: max(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-argmax`, :ref:`label-maximum`, :ref:`label-min`, :ref:`label-sort`.

.. _label-min:

min
---
.. doxygenfunction:: min(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: min(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-argmin`, :ref:`label-minimum`, :ref:`label-max`, :ref:`label-sort`.

.. _label-maximum:

maximum
-------
.. doxygenfunction:: maximum(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: maximum(matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: maximum(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-minimum`, :ref:`label-max`, :ref:`label-min`, :ref:`label-sort`.

.. _label-mean:

mean
----
.. doxygenfunction:: mean(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: mean(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-median`, :ref:`label-variance`, :ref:`label-stddev`.

.. _label-median:

median
------
.. doxygenfunction:: median(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: median(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-variance`, :ref:`label-stddev`.

.. _label-minimum:

minimum
-------
.. doxygenfunction:: minimum(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction:: minimum(matrix<R> const &A, S B)
   :project: castor

.. doxygenfunction:: minimum(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-maximum`, :ref:`label-min`, :ref:`label-max` :ref:`label-sort`.

.. _label-mtimes:

mtimes
------
.. doxygenfunction:: mtimes(R A, matrix<S> const &B)
   :project: castor
.. doxygenfunction::  mtimes(matrix<R> const &A, S B)
   :project: castor
.. doxygenfunction:: mtimes(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-tgemm-naive`, :ref:`label-kron`.

.. _label-norm:

norm
----
.. doxygenfunction:: norm(matrix<S> const &A, std::string typ = "2")
   :project: castor
.. doxygenfunction:: norm(matrix<S> const &A, std::string typ, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-sum`.

.. _label-prod:

prod
----
.. doxygenfunction:: prod(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: prod(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-sum`, :ref:`label-diff`.

.. _label-setdiff:

setdiff
-------
.. doxygenfunction:: setdiff
   :project: castor

See :ref:`label-argsetdiff`, :ref:`label-intersect`, :ref:`label-unique`, :ref:`label-union2`.

.. _label-sort:

sort
----
.. doxygenfunction:: sort
   :project: castor

See :ref:`label-argsort`, :ref:`label-min`, :ref:`label-max`.

.. _label-stddev:

stddev
------
.. doxygenfunction:: stddev(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: stddev(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-mean`, :ref:`label-median`, :ref:`label-variance`.

.. _label-sum:

sum
---
.. doxygenfunction:: sum(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: sum(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-prod`, :ref:`label-diff`, :ref:`label-cumsum`.


.. _label-tgemm-naive:

tgemm
-----
.. doxygenfunction:: tgemm(P alpha, matrix<Q> const &A, matrix<R> const &B, S beta, matrix<T> &C)
   :project: castor


See :ref:`label-mtimes`, :ref:`label-tgemm-blas`.

.. _label-union2:

union2
------
.. doxygenfunction:: union2
   :project: castor

See :ref:`label-intersect`, :ref:`label-setdiff`, :ref:`label-unique`.

.. _label-unique:

unique
------
.. doxygenfunction:: unique
   :project: castor

See :ref:`label-argunique`, :ref:`label-intersect`, :ref:`label-setdiff`, :ref:`label-union2`.

.. _label-variance:

variance
--------
.. doxygenfunction:: variance(matrix<T> const &A)
   :project: castor
.. doxygenfunction:: variance(matrix<T> const &A, int dim)
   :project: castor

See :ref:`label-max`, :ref:`label-min`, :ref:`label-median`, :ref:`label-mean`, :ref:`label-stddev`.
