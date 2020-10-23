.. _label-tools:

Tools
+++++

In this section, we describe some useful tools provided with the **castor** framework. For example, it is possible to produce formatted output with :ref:`label-disp`. It is also possible to get :ref:`label-help` or measure the execution time with :ref:`label-tic` and :ref:`label-toc`.


.. _label-disp:

disp
----
.. doxygenfunction:: disp(matrix<T> const &A, int info = 2, std::ostream &flux = std::cout, std::size_t m = 3, std::size_t n = 3)
   :project: castor

.. doxygenfunction:: disp(T A, int info = 1, std::ostream &flux = std::cout)
   :project: castor

See :ref:`label-help`, :ref:`label-error`, :ref:`label-warning`.

.. _label-help:

help
-----
.. doxygenfunction:: help
   :project: castor

See :ref:`label-disp`, :ref:`label-error`, :ref:`label-warning`.

.. _label-error:

error
-----
.. doxygenfunction:: error
   :project: castor

See :ref:`label-warning`, :ref:`label-disp`, :ref:`label-help`.

.. _label-tic:

tic
---
.. doxygenfunction:: tic
   :project: castor

See :ref:`label-toc`.

.. _label-toc:

toc
---
.. doxygenfunction:: toc
   :project: castor

See :ref:`label-tic`.

.. _label-warning:

warning
-------
.. doxygenfunction:: warning
   :project: castor

See :ref:`label-error`, :ref:`label-disp`, :ref:`label-help`.
