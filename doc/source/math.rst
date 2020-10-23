.. _label-mathematical-functions:

Mathematical functions
++++++++++++++++++++++

The functions described below implement common mathematical functions which can be applied element-wise to the elements of a ``matrix`` such as :ref:`label-cos`, :ref:`label-sqrt`, etc.


.. _label-abs:

abs
---
.. doxygenfunction:: abs
   :project: castor

See :ref:`label-angle`, :ref:`label-sign`.


.. _label-acos:

acos
----
.. doxygenfunction:: acos
   :project: castor

See :ref:`label-acosd`, :ref:`label-cos`.


.. _label-acosd:

acosd
-----
.. doxygenfunction:: acosd
   :project: castor

See :ref:`label-acos`, :ref:`label-cosd`.


.. _label-acosh:

acosh
-----
.. doxygenfunction:: acosh
   :project: castor

See :ref:`label-cosh`.


.. _label-angle:

angle
-----
.. doxygenfunction:: angle
   :project: castor
   
See :ref:`label-abs`.


.. _label-asin:

asin
----
.. doxygenfunction:: asin
   :project: castor

See :ref:`label-sin`, :ref:`label-asind`.


.. _label-asind:

asind
-----
.. doxygenfunction:: asind
   :project: castor

See :ref:`label-sind`, :ref:`label-asin`.


.. _label-asinh:

asinh
-----
.. doxygenfunction:: asinh
   :project: castor

See :ref:`label-sinh`.


.. _label-atan:

atan
----
.. doxygenfunction:: atan
   :project: castor

See :ref:`label-tan`, :ref:`label-atand`.


.. _label-atand:

atand
-----
.. doxygenfunction:: atand
   :project: castor

See :ref:`label-tand`, :ref:`label-atan`.


.. _label-atanh:

atanh
-----
.. doxygenfunction:: atanh
   :project: castor

See :ref:`label-tanh`.


.. _label-ceil:

ceil
----
.. doxygenfunction:: ceil
   :project: castor

See :ref:`label-floor`, :ref:`label-round`.


.. _label-conj:

conj
----
.. doxygenfunction:: conj(matrix<float> const &A)
   :project: castor
.. doxygenfunction:: conj(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: conj(matrix<S> const &X)
   :project: castor

See :ref:`label-real`, :ref:`label-imag`.


.. _label-cos:

cos
---
.. doxygenfunction:: cos
   :project: castor

See :ref:`label-acos`, :ref:`label-cosd`.


.. _label-cosd:

cosd
----
.. doxygenfunction:: cosd
   :project: castor

See :ref:`label-acosd`, :ref:`label-cos`.


.. _label-cosh:

cosh
----
.. doxygenfunction:: cosh
   :project: castor

See :ref:`label-acosh`.


.. _label-deg2rad:

deg2rad
-------
.. doxygenfunction:: deg2rad(T x)
   :project: castor
.. doxygenfunction:: deg2rad(matrix<T> const &X)
   :project: castor

See :ref:`label-rad2deg`.


.. _label-exp:

exp
---
.. doxygenfunction:: exp
   :project: castor

See :ref:`label-log`, :ref:`label-log10`.


.. _label-floor:

floor
-----
.. doxygenfunction:: floor
   :project: castor

See :ref:`label-ceil`, :ref:`label-round`.


.. _label-imag:

imag
----
.. doxygenfunction:: imag(matrix<float> const &A)
   :project: castor
.. doxygenfunction:: imag(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: imag(matrix<S> const &X)
   :project: castor

See :ref:`label-real`, :ref:`label-conj`, :ref:`label-angle`, :ref:`label-abs`.


.. _label-log:

log
---
.. doxygenfunction:: log
   :project: castor

See :ref:`label-log2`, :ref:`label-log10`, :ref:`label-exp`.


.. _label-log2:

log2
----
.. doxygenfunction:: log2
   :project: castor

See :ref:`label-log`, :ref:`label-log10`, :ref:`label-exp`.


.. _label-log10:

log10
-----
.. doxygenfunction:: log10
   :project: castor

See :ref:`label-log`, :ref:`label-log2`, :ref:`label-exp`.


.. _label-pow:

pow
---
.. doxygenfunction:: pow(R x, matrix<S> const &Y)
   :project: castor
.. doxygenfunction:: pow(matrix<R> const &X, S y)
   :project: castor
.. doxygenfunction:: pow(matrix<R> const &X, matrix<S> const &Y)
   :project: castor

See :ref:`label-exp`, :ref:`label-log`.


.. _label-rad2deg:

rad2deg
-------
.. doxygenfunction:: rad2deg(T x)
   :project: castor
.. doxygenfunction:: rad2deg(matrix<T> const &X)
   :project: castor

See :ref:`label-deg2rad`.


.. _label-real:

real
----
.. doxygenfunction:: real(matrix<float> const &A)
   :project: castor
.. doxygenfunction:: real(matrix<double> const &A)
   :project: castor
.. doxygenfunction:: real(matrix<S> const &X)
   :project: castor

See :ref:`label-imag`, :ref:`label-conj`, :ref:`label-angle`, :ref:`label-abs`.


.. _label-round:

round
-----
.. doxygenfunction:: round
   :project: castor

See :ref:`label-floor`, :ref:`label-ceil`.


.. _label-sign:

sign
----
.. doxygenfunction:: sign
   :project: castor

See :ref:`label-abs`.


.. _label-sin:

sin
---
.. doxygenfunction:: sin
   :project: castor

See :ref:`label-asin`, :ref:`label-sind`.


.. _label-sind:

sind
----
.. doxygenfunction:: sind
   :project: castor

See :ref:`label-asind`, :ref:`label-sin`.


.. _label-sinh:

sinh
----
.. doxygenfunction:: sinh
   :project: castor

See :ref:`label-asinh`.


.. _label-sqrt:

sqrt
----
.. doxygenfunction:: sqrt
   :project: castor

See :ref:`label-pow`.


.. _label-tan:

tan
---
.. doxygenfunction:: tan
   :project: castor

See :ref:`label-atan`, :ref:`label-tand`.


.. _label-tand:

tand
----
.. doxygenfunction:: tand
   :project: castor

See :ref:`label-atand`, :ref:`label-tan`.


.. _label-tanh:

tanh
----
.. doxygenfunction:: tanh
   :project: castor

See :ref:`label-atanh`.