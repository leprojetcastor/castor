Replicator mutator
==================

*Shared by Antoine Rideau thanks to Gael Raoul*

On this page you will find how to simulate using **Castor** the replicator mutator equation

.. math:: 

    \partial_{t}u = \underbrace{\sigma \Delta_{x}u}_{mutations} + \underbrace{(x - \bar{x}(t))u}_{replication} \text{ , } t > 0, 

| where :
|    :math:`x \in \mathbb{R}` : a one dimension fitness space,
|    :math:`u(t,x)` : density of a population at time :math:`t` and per unit of fitness,
|    :math:`\bar{x}(t):= \int_{\mathbb{R}}xu(x,t)dx` : mean fitness at time :math:`t` .


We consider a constant population, so

.. math::

    \int_{\mathbb{R}}u(x,t)dx = 1 .



Numeric simulation
------------------

|   Our population gathers ``N`` individuals, each characterized by their fitness :math:`x_{i}`.
|   Ths population will evolve during ``gmax`` generations separated by ``dt`` .

.. code-block:: c++

    // Parameters
    int N = pow(10, 3); // Population
    int gmax = 5000;    // Number of generations
    double dt = 0.01;   // Time discretization
    double sigma = 0.5; // Mutation


Initially, the population is distributed following a Gaussian using ``randn`` .

.. code-block:: c++

    // Initial data
    auto population = randn(1, N);

See :ref:`label-randn` . 

Each generation :

#. Each individuals has a probability :math:`\mathbb{P} = (x_{i})_{+} \times \Delta t` ,where :math:`(x_{i})_{+}` stands for the positive part of :math:`x_{i}` , to give birth to a child who will inherit a fitness of :math:`x_{i} + X` with :math:`X \sim \mathcal{N}(0, \sigma^2 \Delta t)` .

.. code-block:: c++

    // Reproduction
    auto P = dt * max(cat(1, zeros(size(population)), population), 1); 
    auto R = rand(size(population));
    for (int i = 0; i < N; i++)
    {
        if (R(i) < P(i))
        {
            population = cat(2, population, population(i) + sigma * sqrt(dt) * normalDistribution(N * std::rand() / tmp));
        }
    }

See :ref:`label-cat` , :ref:`label-size` , :ref:`label-rand` . 

#. We uniformly choose ``N`` individuals to survive.

.. code-block:: c++

    // Selection
    matrix<std::size_t> killIndex, tmp;
    std::tie(killIndex, tmp) = argunique<std::size_t>(size(population, 2) * rand(1, size(population, 2)));
    killIndex = eval(killIndex(range(0, size(population, 2) - N))); // Index of the individuals who will be killed 
    matrix<std::size_t> Index = setdiff(range(0, size(population, 2)), killIndex); // Index of surviving individuals by taking the complementary of killIndex
    population = eval(population(Index));

See :ref:`label-argunique` , :ref:`label-setdiff` .

Code
----

.. code-block:: c++

    #include "castor/matrix.hpp"
    #include "castor/graphics.hpp"

    using namespace castor;

    int main(int argc, char const *argv[])
    {
        // Parameters
        int N = pow(10, 3); // Population
        int gmax = 10000;   // Number of generations
        double dt = 0.01;   // Time disretization
        double sigma = 0.5; // Mutation
        double tmp = RAND_MAX;

        // Initial data
        auto population = randn(1, N);
        auto normalDistribution = randn(1, N);

        figure fig;

        tic();
        for (int g = 1; g <= gmax; g++)
        {

            std::cout << "---------- Generation " << g << " ----------" << endl;

            // Reproduction
            auto P = dt * max(cat(1, zeros(size(population)), population), 1); 
            auto R = rand(size(population));
            for (int i = 0; i < N; i++)
            {
                if (R(i) < P(i))
                {
                    population = cat(2, population, population(i) + sigma * sqrt(dt) * eval(rand(1, N)(N * rand(1)))); // On a pas le même resultat
                }
            }

            // Selection
            matrix<std::size_t> killIndex, tmp;
            std::tie(killIndex, tmp) = argunique<std::size_t>(size(population, 2) * rand(1, size(population, 2)));
            killIndex = eval(killIndex(range(0, size(population, 2) - N))); // Index of the individuals who will be killed 
            matrix<std::size_t> Index = setdiff(range(0, size(population, 2)), killIndex); // Index of surviving individuals by taking the complementary of killIndex
            population = eval(population(Index));
            
            // Plot
            plot(fig, population, g * dt * ones(size(population)), {"b"});
        }
        toc();

        drawnow(fig);

        return 0;
    }


Reference
---------

https://www.cirm-math.fr/RepRenc/1315/PDFfiles1315.pdf







