Boids
=====

*Shared by Antoine Rideau*

| On this page you will find how to simulate using **Castor** the flocking behaviour of birds.
| Boids - contraction of "bird-oid object" - refers to an artificial life program, developed by Craug Reynolds in 1986, simulating the flocking behaviour of birds.
|
| The complex flocking behaviours arise from three simple rules describing the interaction between each boid:
|       **Separation:** steer to avoid crowding local flockmates.
|       **Alignement:** steer towards the average heading of local flockmates.
|       **Cohesion:** steer to move toward the average position of local flockmates.




The flock is composed of ``N`` boids and will fly during ``nt`` steps within a perimeter described by ``dimension``

.. code-block:: c++

    // Parameters
    int N = 100;  // Number of boids
    int nt = 300; // Number of time steps
    matrix<> dimension = {50, 50}; // Permieter lengths


| Each boids in the flock is characterized by its position :math:`(x,y)` and its velocity vector :math:`\overrightarrow{v}=(dx,dy)` . Each parameter is stored in a ``matrix<>`` and then gathered inside the ``std::vector`` ``Boids``.
| Furthermore the matrix ``Out`` contains these parameters at each steps in order to output them later on for a satisfaying visualization.

.. code-block:: c++

    // Initialization
    std::vector<matrix<>> Boids = { rand(1, N,true) * dimension(0),             // Boids x position
                                    rand(1, N) * dimension(1),                  // Boids y position
                                    rand(1, N,true) * maxSpeed - maxSpeed / 2,  // Boids x velocity : dx
                                    rand(1, N) * maxSpeed - maxSpeed / 2};      // Boids y velocity : dy
    matrix<> Out = cat(1, Boids[0], cat(1, Boids[1], cat(1, Boids[2], Boids[3])));

See :ref:`label-rand` , :ref:`label-cat` .


.. figure:: img/boids.png
    :width: 200
    :align: center
    :figclass: align-center
    
    Neighbourhood of a boid (blue) : zone of separation (red), zone of alignement (yellow) and zone of cohesion (green).

So as to apply the set of rules to each boid, it is needed to determine this boid's flockmates depending on the distance between this boid and the others boids.
This is done by the function ``boidsWithinDistance`` which returns the index of the boids who are inside the circle of ``distance`` radius  around the boid of index ``boidIndex``.

.. code-block:: c++

    matrix<size_t> boidsWithinDistance(std::vector<matrix<>> &Boids, int boidIndex, double distance)
    {
        matrix<size_t> I(boidIndex);
        return setdiff(find(sqrt((Boids[0](boidIndex) - Boids[0]) * (Boids[0](boidIndex) - Boids[0]) + (Boids[1](boidIndex) - Boids[1]) * (Boids[1](boidIndex) - Boids[1])) < distance), I);
    }

See :ref:`label-setdiff` , :ref:`label-find` .


Behaviour
---------

Borders
^^^^^^^

| There are two boundaries conditions that can be applied :

1. Boids must stay within the broders.

    .. code-block:: c++

        // Constrain boids to within the borders. If it gets too close to an edge,nudge it back in and reverse its direction.
        void stayWithinBorders(std::vector<matrix<>> &Boids, matrix<> border, double margin, double turnFactor)
        {
            auto xTooLow = find(Boids[0] < margin);
            auto xTooHigh = find(Boids[0] > border(0) - margin);
            auto yTooLow = find(Boids[1] < margin);
            auto yTooHigh = find(Boids[1] > border(1) - margin);

            if (numel(xTooLow) > 0){Boids[2](xTooLow) = eval(Boids[2](xTooLow)) + turnFactor;}
            if (numel(xTooHigh) > 0){Boids[2](xTooHigh) = eval(Boids[2](xTooHigh)) - turnFactor;}
            if (numel(yTooLow) > 0){Boids[3](yTooLow) = eval(Boids[3](yTooLow)) + turnFactor;}
            if (numel(yTooHigh > 0)){Boids[3](yTooHigh) = eval(Boids[3](yTooHigh)) - turnFactor;}
        }

See :ref:`label-find` , :ref:`label-numel` , :ref:`label-view`  .

2. Like Pacman, boids who go out by one side are wrap around to the other side.

    .. code-block:: c++

        //  Checks if boids go out of the window and if so, wraps them around to the other side.
        void wrapBorders(std::vector<matrix<>> &Boids, matrix<> border)
        {
            auto xTooLow = find(Boids[0] < 0);
            auto xTooHigh = find(Boids[0] > border(0));
            auto yTooLow = find(Boids[1] < 0);
            auto yTooHigh = find(Boids[1] > border(1));

            if (numel(xTooLow) > 0){Boids[0](xTooLow) = eval(Boids[0](xTooLow)) + border(0);}
            if (numel(xTooHigh) > 0){Boids[0](xTooHigh) = eval(Boids[0](xTooHigh)) - border(0);}
            if (numel(yTooLow) > 0){Boids[1](yTooLow) = eval(Boids[1](yTooLow)) + border(1);}
            if (numel(yTooHigh > 0)){Boids[1](yTooHigh) = eval(Boids[1](yTooHigh)) - border(1);}
        }

See :ref:`label-find` , :ref:`label-numel` , :ref:`label-view`  .


Separation
^^^^^^^^^^

For each boid, once flockmates within ``separationDistance`` are known, velocity vector is ajusted toward the opposite directions of flockmates' positions. 

.. code-block:: c++

    // Boids try to keep a small distance away from other boids
    void separation(std::vector<matrix<>> &Boids, double separationDistance, double separationFactor)
    {
        for (int i = 0; i < numel(Boids[0]); i++)
        {
            auto FlockMates = boidsWithinDistance(Boids, i, separationDistance);
            if (numel(FlockMates) > 0)
            {

                Boids[2](i) += sum(Boids[0](i) - eval(Boids[0](FlockMates))) * separationFactor;
                Boids[3](i) += sum(Boids[1](i) - eval(Boids[1](FlockMates))) * separationFactor;
            }
        }
    }

See :ref:`label-sum` , :ref:`label-view`  .



Alignement
^^^^^^^^^^

For each boid, once flockmates within ``alignementDistance`` are known, velocity vector is ajusted toward the average direction of flockmates' velocity vector. 

.. code-block:: c++

    // Boids try to match velocity with near boids
    void alignement(std::vector<matrix<>> &Boids, double alignementDistance, double alignementFactor)
    {
        for (int i = 0; i < numel(Boids[0]); i++)
        {
            auto FlockMates = boidsWithinDistance(Boids, i, alignementDistance);
            if (numel(FlockMates) > 0)
            {
                Boids[2](i) += (sum(eval(Boids[2](FlockMates))) / numel(FlockMates)) * alignementFactor;
                Boids[3](i) += (sum(eval(Boids[3](FlockMates))) / numel(FlockMates)) * alignementFactor;
            }
        }
    }

See :ref:`label-sum` , :ref:`label-numel` , :ref:`label-view`  .

Cohesion
^^^^^^^^

For each boid, once flockmates within ``cohesionDistance`` are known, velocity vector is ajusted toward the average position of flockmates.

.. code-block:: c++

    // Boids try to fly towards the centre of mass of neighbouring boids
    void cohesion(std::vector<matrix<>> &Boids, double cohesionDistance, double cohesionFactor)
    {    
        for (int i = 0; i < numel(Boids[0]); i++)
        {
            auto FlockMates = boidsWithinDistance(Boids, i, cohesionDistance);
            if (numel(FlockMates) > 0)
            {
                Boids[2](i) += ((sum(eval(Boids[0](FlockMates))) / numel(FlockMates)) - Boids[0](i)) * cohesionFactor;
                Boids[3](i) += ((sum(eval(Boids[1](FlockMates))) / numel(FlockMates)) - Boids[1](i)) * cohesionFactor;
            }
        }
    }

See :ref:`label-sum` , :ref:`label-numel` , :ref:`label-view`  .

Speed limitation
^^^^^^^^^^^^^^^^

Once the rules are applied on the boids, their speed is limited to ``maxSpeed`` .

.. code-block:: c++

    void limitSpeed(std::vector<matrix<>> &Boids, double maxSpeed)
    {
        auto Speed = sqrt(Boids[2] * Boids[2] + Boids[3] * Boids[3]);
        matrix<> I = find(Speed > maxSpeed); //Index where velocity > maxSpeed

        if (numel(I) > 0)
        {
            Boids[2](I) = (eval(Boids[2](I)) / eval(Speed(I))) * maxSpeed;
            Boids[3](I) = (eval(Boids[3](I)) / eval(Speed(I))) * maxSpeed;
        }
    }

See :ref:`label-sqrt` , :ref:`label-find` , :ref:`label-view`  .

Visualization
-------------

| In order to have a satisfaying animation of the flying boids Python wille be used.
| Beforehand, the data stored in ``Out`` are written in text file.

.. code-block:: c++

    writetxt("../", "data.txt", Out);


.. code-block:: python

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import matplotlib.patches as patches
    import numpy as np

    data = np.loadtxt("./data.txt", skiprows=1)
    # print(data)

    nt = int(data.shape[0]/4)
    N = int(data.shape[1])

    xs = [data[i*4, :] for i in range(nt)]
    ys = [data[i*4+1, :] for i in range(nt)]
    dxs = [data[i*4+2, :] for i in range(nt)]
    dys = [data[i*4+3, :] for i in range(nt)]

    # Visu initialization
    fig, ax = plt.subplots()
    p = patches.Rectangle((0, 0), 50, 50, fill=True)
    ax.add_artist(p)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 50)
    ax.set_aspect('equal')

    # Length of velocity vector to then normalize its representation
    D = [np.sqrt((dxs[0][i])**2+(dys[0][i])**2) for i in range(N)]
    boids = [ax.annotate("", xy=(xs[0][i]+dxs[0][i]/D[i], ys[0][i] + dys[0][i]/D[i]), xytext=(xs[0][i], ys[0][i]), arrowprops={
                         "facecolor": "red", 'arrowstyle': 'wedge'})for i in range(N)]


    def animate(frame):
        for i in range(N):
            d = np.sqrt((dxs[frame][i])**2+(dys[frame][i])**2)
            pos = np.array([xs[frame][i], ys[frame][i]])

            boids[i].set_position(pos)
            boids[i].xy = pos + (dxs[frame][i]/d, dys[frame][i]/d)

        return boids


    # Creating the Animation object
    ani = animation.FuncAnimation(
        fig, animate, nt, interval=40, blit=True, repeat_delay=1000)
    plt.show()

.. raw:: html

    <video autoplay controls width="100%">

    <source src="./_static/boidsflight.mp4"
            type="video/mp4">

    Sorry, your browser doesn't support embedded videos.
    </video>

|                           Flight of 100 boids in a 50-by-50 square.

Code
----

Here is all the code at once, without the functions written above :


.. code-block:: c++

    int main(int argc, char const *argv[])
    {
        // Parameters
        int N = 100;  // Number of boids
        int nt = 300; // Number of time steps
        matrix<> dimension = {50, 50};

        double margin = 2.;
        double turnFactor = 1.;
        double maxSpeed = 1.5;
        double separationDistance = 2.0;
        double separationFactor = 0.05;
        double alignementDistance = 7.5;
        double alignementFactor = 0.10;
        double cohesionDistance = 8.5;
        double cohesionFactor = 0.03;

        // Initialization
        std::vector<matrix<>> Boids = { rand(1, N,true) * dimension(0),              // Boids x position
                                        rand(1, N) * dimension(1),              // Boids y position
                                        rand(1, N,true) * maxSpeed - maxSpeed / 2,   // Boids x velocity : dx
                                        rand(1, N) * maxSpeed - maxSpeed / 2};  // Boids y velocity : dy
        matrix<> Out = cat(1, Boids[0], cat(1, Boids[1], cat(1, Boids[2], Boids[3])));

        disp(Out);

        tic();
        for (int t = 0; t < nt; t++)
        {
            separation(Boids, separationDistance, separationFactor);
            alignement(Boids, alignementDistance, alignementFactor);
            cohesion(Boids, cohesionDistance, cohesionFactor);
            limitSpeed(Boids, maxSpeed);
            stayWithinBorders(Boids, dimension, margin, turnFactor);

            // Update positions
            Boids[0] += Boids[2];
            Boids[1] += Boids[3];

            // wrapBorders(Boids, dimension);

            for (int i = 0; i < 4; i++)
            {
                Out = cat(1, Out, Boids[i]);
            }
        }
        toc();
        // disp(Out);

        writetxt("../", "data.txt", Out);

        return 0;
    }



References  
----------

https://eater.net/boids

https://en.wikipedia.org/wiki/Boids

http://www.red3d.com/cwr/boids/

http://www.vergenet.net/~conrad/boids/pseudocode.html