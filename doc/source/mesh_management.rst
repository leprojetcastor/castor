.. _label-mesh-management:

Mesh management
+++++++++++++++

Functions helpers to create simple volume or surface meshes from sets of nodes.

.. _label-tetboundary:

tetboundary
-----------
.. doxygenfunction:: tetboundary(matrix<std::size_t> const &tet, matrix<T> const &vtx)
    :project: castor

See :ref:`label-tetdelaunay`.


.. _label-tetdelaunay:

tetdelaunay
-----------
.. doxygenfunction:: tetdelaunay(matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z)
    :project: castor

See :ref:`label-tetboundary`.


.. _label-tridelaunay:

tridelaunay
-----------
.. doxygenfunction:: tridelaunay(matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z = {})
    :project: castor
