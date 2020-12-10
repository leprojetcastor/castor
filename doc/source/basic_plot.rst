.. _label-basic-plot:


Basic plot
++++++++++

These functions allow to display curves or values of a ``matrix``.



.. _label-imagesc:

imagesc
-------
.. doxygenfunction:: imagesc(figure &fig, matrix<T> const &M)
    :project: castor

.. _label-plot:

plot
----
.. doxygenfunction:: plot(figure &fig, matrix<T> const &X, matrix<T> const &Y, std::vector<std::string> const &style = {""}, std::vector<std::string> const &label = {""})
    :project: castor

See :ref:`label-plot3`

.. _label-plot3:

plot3
-----
.. doxygenfunction:: plot3(figure &fig, matrix<T> const &X, matrix<T> const &Y, matrix<T> const &Z, std::string const &style = "")
    :project: castor

See :ref:`label-plot`

.. _label-spy:

spy
---
.. doxygenfunction:: spy
    :project: castor