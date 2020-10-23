.. _label-builders:

Builders
++++++++

In this section, the user will find all the possible builders other than the constructors for the ``matrix<>`` object. For the latter, please refer to constructors list in the :ref:`class matrix description <label-class-matrix>`. For example, :ref:`label-eye` will return the idendity matrix while :ref:`label-ones` will return a matrix filled with ones.

.. _label-colon:

colon
-----
.. doxygenfunction:: colon(U j, V k)
   :project: castor
.. doxygenfunction:: colon(U j, V i, W k)
   :project: castor

See :ref:`label-linspace`.

.. _label-diag:

diag
----
.. doxygenfunction:: diag
   :project: castor

See :ref:`label-eye`.

.. _label-eye:

eye
---
.. doxygenfunction:: eye(matrix<std::size_t> const &S)
   :project: castor
.. doxygenfunction:: eye(std::size_t m, long n = -1)
   :project: castor

See :ref:`label-ones`, :ref:`label-zeros`, :ref:`label-rand`.

.. _label-linspace:

linspace
--------
.. doxygenfunction:: linspace
   :project: castor

See :ref:`label-logspace`, :ref:`label-colon`.

.. _label-logspace:

logspace
--------
.. doxygenfunction:: logspace
   :project: castor

See :ref:`label-linspace`, :ref:`label-colon`, :ref:`label-log`.

.. _label-rand:

rand
----
.. doxygenfunction:: rand(matrix<std::size_t> const &S, bool seed = false)
   :project: castor
.. doxygenfunction:: rand(std::size_t m, long n = -1, bool seed = false)
   :project: castor

See :ref:`label-zeros`, :ref:`label-eye`, :ref:`label-ones`.

.. _label-ones:

ones
----
.. doxygenfunction:: ones(matrix<std::size_t> const &S)
   :project: castor
.. doxygenfunction:: ones(std::size_t m, long n = -1)
   :project: castor

See :ref:`label-zeros`, :ref:`label-eye`, :ref:`label-rand`.

.. _label-zeros:

zeros
-----
.. doxygenfunction:: zeros(matrix<std::size_t> const &S)
   :project: castor
.. doxygenfunction:: zeros(std::size_t m, long n = -1)
   :project: castor

See :ref:`label-ones`, :ref:`label-eye`, :ref:`label-rand`.
