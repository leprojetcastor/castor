
.. _label-class-matrix:

Class matrix
+++++++++++++

The **castor** framework implements its own templatized class ``matrix<T>`` where ``T`` can be *for example* a ``float``, ``int``, ``std::complex<double>``, etc. At its core, it is built over a ``std::vector<T>`` which holds the values. The class ``matrix`` itself provides many useful functions and operators (addition, multiplication, indexing, etc.). It is designed such that it should feel like using `Matlab <https://www.mathworks.com>`_ or `Numpy arrays <htttps://www.numpy.org>`_. The user will find here all the available constructors. Specific builders (:ref:`label-ones`, :ref:`label-eye`, etc.) may be found at :ref:`label-builders`.

.. doxygenclass:: castor::matrix
   :project: castor
   :members: 

See :ref:`label-view`, :ref:`label-cview`.

