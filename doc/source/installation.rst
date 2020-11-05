.. _label-installation:

Installation
============

Download header files
+++++++++++++++++++++

**(the easy way!)**

The simplest way to get the **castor** library is to download the `latest version of header files <https://gitlab.labos.polytechnique.fr/leprojetcastor/castor/-/raw/master/include/matrix?inline=false>`_ and include them during the compilation of your c++ program using the library, see :ref:`label-compilation` in the :ref:`label-basic` section.

From git repository with CMake
++++++++++++++++++++++++++++++

**(the not-too-complicated way)**

You can also install the **castor** library from source with ``CMake``.

On Linux and macOS platforms :

.. code::

    $ git clone git@gitlab.labos.polytechnique.fr:leprojetcastor/castor.git 
    $ cd castor
    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_INSTALL_PREFIX=path/to/install/directory ..
    $ make install

``path/to/install/directory`` is the absolute path to the folder where **castor** header files are installed. ``CMake`` assumes this folder contains ``include`` subfolder. By default, ``CMAKE_INSTALL_PREFIX`` is set to ``/usr/local`` and header files are installed in ``/usr/local/include/castor`` directory. 

The linear algebra part and the visualization part of the library depend respectively on an optimized BLAS library (like OpenBLAS and MKL) and the VTK library. ``CMake`` tries to detect these required libraries in your system, if ``Cmake`` can not find these dependencies, you can indicate the paths to them by using ``-DCMAKE_PREFIX_PATH``:

.. code::

    $ cmake -DCMAKE_INSTALL_PREFIX=path/to/install/directory -DCMAKE_PREFIX_PATH="/path/to/optimized/blas;/path/to/vtk/" ..   

Installing the dependencies
+++++++++++++++++++++++++++

Installing BLAS and LAPACK
..........................

BLAS and LAPACK may be obtained through multiple channels:

 - the first possibility it to download the binaries for `OpenBLAS <https://www.openblas.net/>`_. On Ubuntu, they can be obtained with 

    .. code:: text

        sudo apt install libopenblas-dev

    and on MacOS with

    .. code::

        brew install openblas

 - the second possibility is to download the ``Intel MKL`` framework which is proprietary but freely available for non-commercial use. Binaries can be obtained directly from the Intel website (but it will require the creation of an account) or using `Anaconda <https://www.anaconda.com/>`_. One of the advantages of the ``MKL`` is that it is deeply optimized for those using Intel CPUs. For AMD CPU users, it is a bit trickier and we refer to `https://danieldk.eu/Posts/2020-08-31-MKL-Zen.html <https://danieldk.eu/Posts/2020-08-31-MKL-Zen.html>`_.


**Remark:** on MacOS, the ``vecLib`` framework shipped with the Apple Command Line Tools is an optimized implementation of BLAS and LAPACK. ``vecLib`` will be detected automatically by ``cmake`` if no other BLAS/LAPACK implementations are found.


Installing VTK
..............

The installation of VTK can be more complicated, especially since ``VTK 8.x.x`` is required. 

On MacOS, it should be installed using Homebrew with the command 

.. code:: text

    brew install vtk@8.2

Please note that ``brew install vtk`` (assuming ``Homebrew`` is up-to-date) will install a version from the ``9.x.x`` branch (or later) which is currently not compatible. Moreover, since the path to ``VTK 8.2`` will not be found automatically by ``cmake``, it should be given using the ``-DCMAKE_PREFIX_PATH=/path/to/vtk@8.2/`` flag as explained at the beginning of this page.

On Ubuntu, VTK 8.2 must be built from source as the version available in the repositories is currently ``VTK 7.x.x``. The process is described below:

 - first, download the source code at `https://vtk.org/ <https://vtk.org/>`_ and extract the archive.

 - open a terminal and install ``freeglut3-dev``.

    .. code:: text

        sudo apt install freeglut3-dev

 - go to the main directory of VTK, create a ``build`` directory and go to the newly created directory.

    .. code:: text

        cd /path/to/VTK/main
        mkdir build && cd build


 - call the following ``cmake`` command.

    .. code:: text

        cmake -DCMAKE_BUILD_TYPE=Release ..

    If you wish to change the default install directory, add the following flag to the command above: ``-DCMAKE_INSTALL_PREFIX=/your/path/``. In order to use VTK, you will need to give ``/your/path/`` when configuring your project with ``cmake -DCMAKE-DCMAKE_PREFIX_PATH=/your/path/ ..``, see above.

    **Remark:** if you are not familiar with this process, do not modify the default installation folder.

 - compile and install.

    .. code:: text

        make
        sudo make install

    **Remark:** if you have a CPU with ``N`` cores (not *threads*), you can accelerate the compilation of VTK with

    .. code:: text

        make -jN
        sudo make install

    **Remark n**:math:`^{o} 2` **:** the compilation (the ``make -jN`` command) will take some time, so you can go grab yourself a cup of tea or coffee...

The binaries and the headers are, normally, placed respectively in the ``/usr/local/lib/`` and ``/usr/local/include`` folders and should be found automatically by ``cmake``.

You can also create a file ``install_vtk.sh`` (or whatever the name you wish, but with the ``.sh`` extension) with the following content

.. code:: text

    sudo apt install freeglut3-dev
    wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
    cd VTK-8.2.0/
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/your/vtk/install/folder -DCMAKE_BUILD_TYPE=Release ..
    make -jN
    sudo make install

where the ``-DCMAKE_INSTALL_PREFIX=...`` option may be omitted if you want to use the default installation folder and ``N`` is set to the number of CPU *cores*. Then, simply call

.. code:: text

    bash install_vtk.sh

to start the installation process.
