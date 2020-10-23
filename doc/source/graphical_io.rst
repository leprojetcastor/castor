.. _label-graphical-io:

Graphical input/output
++++++++++++++++++++++

.. _label-triread:

triread
-------
.. doxygenfunction:: triread(std::string const &path, std::string const &name)
    :project: castor

See :ref:`label-triwrite`.


.. _label-triwrite:

triwrite
--------
.. doxygenfunction:: triwrite(std::string const &path, std::string const &name, matrix<std::size_t> const &tri, matrix<T> const &vtx, matrix<T> const &val = {})
    :project: castor

See :ref:`label-triread`.

