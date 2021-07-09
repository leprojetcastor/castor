N-body problem
==============

*Shared by Antoine Rideau*

| On this page you will find how to simulate using **Castor** the n-body problem with 3 celestial bodies : the Sun, Jupiter and Saturn.
|
| In physics, the n-body problem is the problem of predicting the individual motions of a group of celestial objects interacting with each other gravitationally.
|
| First of all,  Newton's law of gravity says that the gravitational force felt on a planet *P* by the Sun *S* is given by

.. math::

    \overrightarrow{F}_{S\rightarrow P} = - \overrightarrow{F}_{P\rightarrow S} = -\frac{Gm_{S}m_{P}}{d^2}\overrightarrow{u} ,

| where :
|    :math:`G` : gravitational constant,
|    :math:`m_{S}`, :math:`m_{P}` : masses of the Sun and the planet,
|    :math:`d` : distance between the Sun and the planet,
|    :math:`\overrightarrow{u}` : monic vector directed from the Sun to the planet.
|
| The system is composed of the three heavier bodies in the solar system : the Sun, Jupiter and Saturn which will be referred as *S*, *Ju* and *Sa*. This system is supposed isolated and only gravitational forces are applied on the three celestial bodies.
| According to Newton's second principle, with :math:`m` the masses, :math:`\overrightarrow{a}` the acceleration and :math:`\overrightarrow{F}` the forces :

.. math::

    \begin{matrix}
    m_{Ju}\overrightarrow{a}_{Ju} = \overrightarrow{F}_{S\rightarrow Ju}  + \overrightarrow{F}_{Sa\rightarrow Ju} ,
    \\ 
    m_{Sa}\overrightarrow{a}_{Sa} = \overrightarrow{F}_{S\rightarrow Sa}  + \overrightarrow{F}_{Ju\rightarrow Sa} ,
    \\ 
    m_{S}\overrightarrow{a}_{S} = \overrightarrow{F}_{Ju\rightarrow S}  + \overrightarrow{F}_{Sa\rightarrow S} .
    \end{matrix}


+------------+------------------------------------------------------+
|   Bodies   |  Masses                                              |
|            |  (relatively to the Sun)                             |
+============+======================================================+
| Sun        | :math:`m_{0}` = 1.00000597682                        |
+------------+------------+-----------------------------------------+
| Jupiter    | :math:`m_{1}` = 0.000954786104043                    |
+------------+------------+-----------------------------------------+
| Saturn     | :math:`m_{2}` = 0.000285583733151                    |
+------------+------------+-----------------------------------------+
| Gravitational constant :math:`G = 2.95912208286 \times 10^{-4}`   |
+-------------------------------------------------------------------+


Those constants are stored within a data strucuture named ``DATA``

.. code-block:: c++

    struct DATA
    {
        double m0;
        double m1;
        double m2;
        double G;
    };

.. code-block:: c++

    // Parameters
    DATA cst;
    cst.m0 = 1.00000597682;    // Sun's mass
    cst.m1 = 9.54786104043e-4; // Jupiter's mass
    cst.m2 = 2.85583733151e-4; // Saturn's mass
    cst.G = 2.95912208286e-4;  //Gravitation's constant


The position of a object over time is given by :math:`q(t)`.


+------------+-----------------------------+
|   Bodies   |    Initial position (AU)    |
+============+=============================+
|            |               0             |
|            +-----------------------------+
|     Sun    |               0             |
|            +-----------------------------+
|            |               0             |
+------------+-----------------------------+
|            |          −3.5023653         |
|            +-----------------------------+
|  Jupiter   |          −3.8169847         |
|            +-----------------------------+
|            |          −1.5507963         |
+------------+-----------------------------+
|            |           9.0755314         |
|            +-----------------------------+
| Saturn     |          −3.0458353         |
|            +-----------------------------+
|            |          −1.6483708         |
+------------+-----------------------------+


where AU stands for astronomical unit and :math:`1 AU = 1.495 978 707 \times 10^{11}` m.

.. code-block:: c++
    
    matrix<> qini = {0, 0, 0,                               // Sun's initial position
                    -3.5023653, -3.8169847, -1.5507963,     // Jupiter's initial position
                    9.0755314, -3.0458353, -1.6483708};     // Saturn's initial position


| Moreover, using the momentum :math:`\overrightarrow{p} = m\overrightarrow{v} = m\overrightarrow{\dot{q}}(t)` instead of the speed is more practical. 
| Indeed, the system's total momentum is preserved over time :

.. math::

    \overrightarrow{p}_{S}(t) + \overrightarrow{p}_{Ju}(t) + \overrightarrow{p}_{Sa}(t) = Constant

+------------+-----------------------------+
|   Bodies   | Initial momentum (AU/day)   |
+============+=============================+
|            |               0             |
|            +-----------------------------+
|     Sun    |               0             |
|            +-----------------------------+
|            |               0             |
+------------+-----------------------------+
|            |           0.00565429        |
|            +-----------------------------+
|  Jupiter   |           0.00565429        |
|            +-----------------------------+
|            |          −0.00190589        |
+------------+-----------------------------+
|            |           0.00168318        |
|            +-----------------------------+
| Saturn     |           0.00483525        |
|            +-----------------------------+
|            |           0.00192462        |
+------------+-----------------------------+


.. code-block:: c++

    matrix<> pini = {0, 0, 0,                                                           // Sun's momentum
                    0.00565429 * cst.m1, -0.00412490 * cst.m1, -0.00190589 * cst.m1,    // Jupiter's momentum
                    0.00168318 * cst.m2, 0.00483525 * cst.m2, 0.00192462 * cst.m2};     // Saturn's momentum


Time is discretized into ``nt`` steps 

.. code-block:: c++

    // Disretization
    int nt = 1501;
    double dt = (tend - tini) / (nt - 1);
    auto T = linspace(tini, tend, nt);


Scheme
------


Code
----


References
----------

| https://interstices.info/les-planetes-tournent-elles-rond/
|
| https://www.f-legrand.fr/scidoc/srcdoc/numerique/symplectic/verlet/verlet-pdf.pdf