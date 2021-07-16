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


Those constants are stored within a data structure named ``DATA``

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

.. math::

    t_{i} = it \times \delta t \text{ for } it = \left [ \! \left [ 0, nt-1 \right ] \! \right ]

.. code-block:: c++

    // Disretization
    int nt = 1501;
    double dt = (tend - tini) / (nt - 1);
    auto T = linspace(tini, tend, nt);

See :ref:`label-linspace`.

Scheme
------

| A symplectic Euler scheme is used in this simulation because it preserved the energy of the system unlike either forward and backward Euler scheme.
| 
| As the system is conservative the Hamiltonian can be separated in a cinetical part :math:`K(p)` and a potential part :math:`V(q)` :

.. math::

    H(q,p) = K(p) + V(q)


With such a separation, Hamilton equation are given by 

.. math::

    \begin{matrix}
    \displaystyle \frac{\mathrm{d} q}{\mathrm{d} t} = + \frac{\mathrm{d} K}{\mathrm{d} p}
    \\
    \\
    \displaystyle \frac{\mathrm{d} p}{\mathrm{d} t} = - \frac{\mathrm{d} V}{\mathrm{d} q}
    \end{matrix}

which result to the symplectic Euler scheme :

.. math::

    \begin{matrix}
    \displaystyle q_{n+1} = q_{n} + \frac{\mathrm{d} K}{\mathrm{d} p}(p_{n})
    \\
    \\
    \displaystyle p_{n+1} = p_{n} - \frac{\mathrm{d} V}{\mathrm{d} q}(q_{n+1})
    \end{matrix}

.. code-block:: c++

    // Scheme
    auto Q = zeros(nt, numel(qini));
    Q(0, col(Q)) = qini;
    auto P = zeros(nt, numel(pini));
    P(0, col(P)) = pini;
    // Symplectic Euler
    for (int it = 0; it < nt - 1; it++)
    {
        matrix<> q_n = eval(Q(it, col(Q)));
        matrix<> p_n = eval(P(it, col(P)));
        Q(it + 1, col(Q)) = q_n + dt * H_p(cst, p_n);
        P(it + 1, col(P)) = p_n - dt * H_q(cst, eval(Q(it + 1, col(Q))));
    }

See :ref:`label-zeros`, :ref:`label-numel`, :ref:`label-col` , :ref:`label-view`.

In the code, :math:`\displaystyle \frac{\mathrm{d} K}{\mathrm{d} p}(p)` is represented by the function ``H_p`` 

.. code-block:: c++

    matrix<> H_p(DATA cst, matrix<> p)
    {
        auto Hp = zeros(1, 9);
        Hp(range(0, 3)) = eval(p(range(0, 3))) / cst.m0;
        Hp(range(3, 6)) = eval(p(range(3, 6))) / cst.m1;
        Hp(range(6, 9)) = eval(p(range(6, 9))) / cst.m2;
        return Hp;
    }

See :ref:`label-zeros`, :ref:`label-range` , :ref:`label-view`.

and :math:`\displaystyle \frac{\mathrm{d} V}{\mathrm{d} q}(q)` by the function ``H_q``

.. code-block:: c++

    matrix<> H_q(DATA cst, matrix<> q)
    {
        double m0 = cst.m0;
        double m1 = cst.m1;
        double m2 = cst.m2;
        double G = cst.G;
        auto q0 = eval(q(range(0, 3)));
        auto q1 = eval(q(range(3, 6)));
        auto q2 = eval(q(range(6, 9)));
        auto Hq = zeros(1, 9);
        Hq(range(0, 3)) = (G * m0 * m1 * ((q0 - q1) / pow(norm(q0 - q1), 3)) + G * m0 * m2 * ((q0 - q2) / pow(norm(q0 - q2), 3)));
        Hq(range(3, 6)) = (G * m1 * m0 * ((q1 - q0) / pow(norm(q1 - q0), 3)) + G * m1 * m2 * ((q1 - q2) / pow(norm(q1 - q2), 3)));
        Hq(range(6, 9)) = (G * m2 * m0 * ((q2 - q0) / pow(norm(q2 - q0), 3)) + G * m2 * m1 * ((q2 - q1) / pow(norm(q2 - q1), 3)));
        return Hq;
    }

See :ref:`label-range`, :ref:`label-view`.


Post processing
---------------

| In order to have a video of the orbiting planets, Python will be used instead of a simple ``plot`` function.
|
| To do so, first the positions in the matrix ``Q`` are stored using ``writetxt`` .

.. code-block:: c++

    // Output
    writetxt("./", "dataJu.txt", cat(2, eval(Q(row(Q), 3)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 4)) - eval(Q(row(Q), 1))));
    writetxt("./", "dataSa.txt", cat(2, eval(Q(row(Q), 6)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 7)) - eval(Q(row(Q), 1))));

See :ref:`label-writetxt`, :ref:`label-row` .

Then the following Python code shows the beautiful animation.

.. code-block:: python

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import numpy as np
    from collections import deque

    # Data input
    dataJu = np.loadtxt("./build/dataJu.txt")
    dataSa = np.loadtxt("./build/dataSa.txt")

    # Parameters extraction
    nt = int(dataJu[0, 0])

    # Data processing
    dataJu = np.delete(dataJu, 0, 0)
    dataSa = np.delete(dataSa, 0, 0)

    # Visu initialization
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-10, 10), ylim=(-10, 10))
    ax.set_aspect('equal')

    line, = ax.plot([], [], 'o', lw=2)
    traceJu, = ax.plot([], [], ',-', lw=1)
    traceSa, = ax.plot([], [], ',-', lw=1)
    historyJu_x, historyJu_y = deque(maxlen=nt), deque(maxlen=nt)
    historySa_x, historySa_y = deque(maxlen=nt), deque(maxlen=nt)


    def animate(i):
        # Get planets' current positions
        thisx = [0, dataJu[i, 0], dataSa[i, 0]]
        thisy = [0, dataJu[i, 1], dataSa[i, 1]]

        # Clear the trace when the animation loops
        if i == 0:
            historyJu_x.clear()
            historyJu_y.clear()
            historySa_x.clear()
            historySa_y.clear()

        # Add the current position to the trace
        historyJu_x.appendleft(thisx[1])
        historyJu_y.appendleft(thisy[1])
        historySa_x.appendleft(thisx[2])
        historySa_y.appendleft(thisy[2])

        line.set_data(thisx, thisy)  # Update planets' positions
        # Update planets' traces
        traceJu.set_data(historyJu_x, historyJu_y)
        traceSa.set_data(historySa_x, historySa_y)
        return line, traceJu, traceSa


    # Creating the Animation object
    ani = animation.FuncAnimation(
        fig, animate, nt, interval=10, blit=True)
    plt.show()




Code
----

Here is all the code at once, without the functions ``H_q`` and ``H_p``  written above :

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/graphics.hpp"
    #include "castor/linalg.hpp"

    using namespace castor;

    struct DATA
    {
        double m0;
        double m1;
        double m2;
        double G;
    };

    int main(int argc, char const *argv[])
    {
        // Parameters
        DATA cst;
        cst.m0 = 1.00000597682;    // Sun's mass
        cst.m1 = 9.54786104043e-4; // Jupiter's mass
        cst.m2 = 2.85583733151e-4; // Saturn's mass
        cst.G = 2.95912208286e-4;  // Gravitation's constant

        matrix<> qini = {0, 0, 0,                                                       // Sun's initial position
                        -3.5023653, -3.8169847, -1.5507963,                             // Jupiter's initial position
                        9.0755314, -3.0458353, -1.6483708};                             // Saturn's initial position
        matrix<> pini = {0, 0, 0,                                                       // Sun's initial momentum
                        0.00565429 * cst.m1, -0.00412490 * cst.m1, -0.00190589 * cst.m1,// Jupiter's initial momentum
                        0.00168318 * cst.m2, 0.00483525 * cst.m2, 0.00192462 * cst.m2}; // Saturn's initial momentum

        double tini = 0.;
        double tend = 12500.;
        
        // Disretization
        int nt = 1501;
        double dt = (tend - tini) / (nt - 1);
        auto T = linspace(tini, tend, nt); 

        // Symplectic Euler scheme
        auto Q = zeros(nt, numel(qini));    // Matrix of positions over time
        Q(0, col(Q)) = qini;                // Initialization 
        auto P = zeros(nt, numel(pini));    // Matrix of momentum over time
        P(0, col(P)) = pini;                // Initialization
        for (int it = 0; it < nt - 1; it++)
        {
            matrix<> q_n = eval(Q(it, col(Q)));
            matrix<> p_n = eval(P(it, col(P)));
            Q(it + 1, col(Q)) = q_n + dt * H_p(cst, p_n);
            P(it + 1, col(P)) = p_n - dt * H_q(cst, eval(Q(it + 1, col(Q))));
        }

        // Output
        writetxt("./", "dataJu.txt", cat(2, eval(Q(row(Q), 3)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 4)) - eval(Q(row(Q), 1))));
        writetxt("./", "dataSa.txt", cat(2, eval(Q(row(Q), 6)) - eval(Q(row(Q), 0)), eval(Q(row(Q), 7)) - eval(Q(row(Q), 1))));
    
        return 0;
    }

    Orbit of Jupiter (orange) and Saturn (green) around the Sun in the center.

.. raw:: html

    <video controls width="100%">

    <source src="./_static/3body.mp4"
            type="video/mp4">

    Sorry, your browser doesn't support embedded videos.
    </video>



References
----------

| https://interstices.info/les-planetes-tournent-elles-rond/
|
| http://www.unige.ch/~hairer/poly/chap3.pdf