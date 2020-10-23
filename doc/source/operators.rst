.. _label-operators:

Operators
+++++++++

This section describes all the conventional ``C++`` operators which have been overloaded in ``matrix``. To the exception of ``<<``, they all work element-wise.


.. _label-operator<<:

operator<<
----------
.. doxygenfunction:: operator<<(std::ostream &flux, matrix<T> const &A)
   :project: castor

See :ref:`label-disp`.


.. _label-operator!:

operator!
---------
.. doxygenfunction:: operator!
   :project: castor

See :ref:`label-operator&&`, :ref:`label-operator||`.


.. _label-operator&&:

operator&&
----------
.. doxygenfunction:: operator&&(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator&&(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator&&(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator||`, :ref:`label-operator!`.


.. _label-operator||:

operator||
----------
.. doxygenfunction:: operator||(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator||(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator||(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator&&`, :ref:`label-operator!`.


.. _label-operator==:

operator==
----------
.. doxygenfunction:: operator==(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator==(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator==(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator!=`, :ref:`label-operator<=`, :ref:`label-operator>=`.


.. _label-operator!=:

operator!=
----------
.. doxygenfunction:: operator!=(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator!=(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator!=(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator==`, :ref:`label-operator<=`, :ref:`label-operator>=`.


.. _label-operator<=:

operator<=
----------
.. doxygenfunction:: operator<=(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator<=(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator<=(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator==`, :ref:`label-operator>=`, :ref:`label-operator<`.


.. _label-operator<:

operator<
---------
.. doxygenfunction:: operator<(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator<(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator<(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator<=`, :ref:`label-operator>`, :ref:`label-operator>=`.


.. _label-operator>=:

operator>=
----------
.. doxygenfunction:: operator>=(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator>=(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator>=(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator>`, :ref:`label-operator<=`, :ref:`label-operator<`.


.. _label-operator>:

operator>
---------
.. doxygenfunction:: operator>(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator>(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator>(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator>=`, :ref:`label-operator<`, :ref:`label-operator<=`.


.. _label-operator+:

operator+
---------
.. doxygenfunction:: operator+(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator+(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator+(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator-`, :ref:`label-operator*`, :ref:`label-operator/`.


.. _label-operator-:

operator-
---------
.. doxygenfunction:: operator-(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator-(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator-(matrix<R> const &A, matrix<S> const &B)
   :project: castor

.. doxygenfunction:: operator-(matrix<S> const &A)
   :project: castor

See :ref:`label-operator+`, :ref:`label-operator*`, :ref:`label-operator/`.


.. _label-operator*:

operator*
---------
.. doxygenfunction:: operator*(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator*(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator*(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator+`, :ref:`label-operator-`, :ref:`label-operator/`.


.. _label-operator/:

operator/
---------
.. doxygenfunction:: operator/(R A, matrix<S> const &B) 
   :project: castor

.. doxygenfunction:: operator/(matrix<R> const &A, S B) 
   :project: castor

.. doxygenfunction:: operator/(matrix<R> const &A, matrix<S> const &B)
   :project: castor

See :ref:`label-operator+`, :ref:`label-operator-`, :ref:`label-operator*`.