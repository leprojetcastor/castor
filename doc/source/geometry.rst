.. _label-geometry:

Geometry
++++++++

The functions from the **Geometry** category are useful to perform transformations between cartesian and polar or spherical coordinates. It is also possible to compute a cartesian 2-dimensional grid with :ref:`label-meshgrid` or to scatter three-dimensional nodes over a sphere with the :ref:`label-sphere` and :ref:`label-sphere2` (Fibonacci sphere) functions.


.. _label-cart2pol:

cart2pol
--------
.. doxygenfunction:: cart2pol(R const &X, S const &Y)
   :project: castor
.. doxygenfunction:: cart2pol(matrix<R> const &X, matrix<S> const &Y)
   :project: castor

See :ref:`label-pol2cart`, :ref:`label-cart2sph`, :ref:`label-sph2cart`.

.. _label-cart2sph:

cart2sph
--------
.. doxygenfunction:: cart2sph(Q const &X, R const &Y, S const &Z)
   :project: castor
.. doxygenfunction:: cart2sph(matrix<Q> const &X, matrix<R> const &Y, matrix<S> const &Z)
   :project: castor

See :ref:`label-sph2cart`, :ref:`label-cart2pol`, :ref:`label-pol2cart`.

.. _label-idx2sph:

idx2sph
--------
.. doxygenfunction:: idx2sph
   :project: castor

See :ref:`label-sph2idx`, :ref:`label-sphere`.

.. _label-meshgrid:

meshgrid
--------
.. doxygenfunction:: meshgrid
   :project: castor

See :ref:`label-kron`.

.. _label-pol2cart:

pol2cart
--------
.. doxygenfunction:: pol2cart(R const &THE, S const &RHO)
   :project: castor
.. doxygenfunction:: pol2cart(matrix<R> const &THE, matrix<S> const &RHO)
   :project: castor

See :ref:`label-cart2pol`, :ref:`label-cart2sph`, :ref:`label-sph2cart`.

.. _label-sph2cart:

sph2cart
--------
.. doxygenfunction:: sph2cart(Q const &THE, R const &PHI, S const &RHO)
   :project: castor
.. doxygenfunction:: sph2cart(matrix<Q> const &THE, matrix<R> const &PHI, matrix<S> const &RHO)
   :project: castor

See :ref:`label-cart2sph`, :ref:`label-cart2pol`, :ref:`label-pol2cart`.

.. _label-sph2idx:

sph2idx
--------
.. doxygenfunction:: sph2idx
   :project: castor

See :ref:`label-idx2sph`, :ref:`label-sphere`.

.. _label-sphere:

sphere
------
.. doxygenfunction:: sphere
   :project: castor

See :ref:`label-idx2sph`, :ref:`label-sph2idx`.

.. _label-sphere2:

sphere2
-------
.. doxygenfunction:: sphere2
   :project: castor

See :ref:`label-sphere`.