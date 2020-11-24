.. _label-installation:

Installation
============

The simplest way to get the **castor** library is to download the `latest version of header files <https://gitlab.labos.polytechnique.fr/leprojetcastor/castor/-/jobs/artifacts/master/download?job=deploy>`_ and include them during the compilation of your c++ program using the library, see :ref:`label-compilation` in the :ref:`label-basic` section.

For a complete installation integrating check dependencies and examples compilation, we describe below the procedure for **MacOS** and **Linux** distributions like Ubuntu. For the installation on **Windows**, see the dedicated section :ref:`label-install-windows`.

From git repository with CMake
++++++++++++++++++++++++++++++

You can install the **castor** library from source with ``CMake``. The source files of the library is available here `<https://gitlab.labos.polytechnique.fr/leprojetcastor/castor>`_.

On Linux and macOS platforms :

.. code::

    $ git clone https://gitlab.labos.polytechnique.fr/leprojetcastor/castor.git
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


.. _label-install-vtk:

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

.. _label-install-windows:

Installing on Windows 10
++++++++++++++++++++++++

There is *a priori* no easy solution on Windows 10. One possibility is to use *only* the Visual Studio tools (freely available for non-commercial use). The blocking point is the compilation of BLAS/LAPACK as it requires a Fortran compiler which is a complicated topic. Consequently, a possibility would be to use the Intel MKL library (also freely available for non-commercial use, but requires a registration). In order to build VTK, one can follow the recommandations `here <https://vtk.org/Wiki/VTK/Building/Windows>`_. The **castor** framework *could* then be installed in a similar fashion as VTK using ``cmake`` or ``cmake-gui``. 

The solution above has not yet been fully tested and we will rather use the `MSYS2 tools <https://www.msys2.org/>`_. MSYS2 will allow the Unix/MacOS user to work with a familiar self-contained environment within Windows. After installation, start a MSYS2 terminal (for a standard installation, the executable file is ``C:\msys64\mingw64.exe``) and update the database with the following commands:

.. code:: text

    pacman -Syu
    pacman -Su

First, install the build tools, ``GCC``, ``git`` and ``cmake``:

.. code:: text

    pacman -S base-devel
    pacman -S mingw-w64-x86_64-gcc
    pacman -S git 
    pacman -S mingw-w64-x86_64-cmake

It may take a *lot of time*. Now, install the dependencies ``OpenBLAS`` and ``VTK 8.2``:

.. code:: text

    pacman -S mingw-w64-x86_64-openblas
    pacman -S mingw-w64-x86_64-vtk

Note that the current version of ``VTK`` is 8.2. If it is not the case, you will need to compile it from source. Fortunately, it happens in the same fashion as for the Ubuntu case, see :ref:`label-install-vtk`. We are now ready to install **castor**. First, clone the repository, and create a ``build directory``:

.. code:: text

    git clone https://gitlab.labos.polytechnique.fr/leprojetcastor/castor.git
    cd castor
    mkdir build && cd build

Now, let us generate the build files. ``VTK`` should normally be found  automatically but it may not be the case for ``OpenBLAS``. The following command should work:

.. code:: text

    cmake -G"MSYS Makefiles" -DBLAS_LIBRARIES="/mingw64/lib/libopenblas.a" -DLAPACK_LIBRARIES="/mingw64/lib/libopenblas.a" -DCBLAS_INCLUDE_DIR="/mingw64/include/OpenBLAS" ..

**Remark:** The ``-G"MSYS Makefiles"`` is mandatory. Otherwise, ``cmake`` could try to generate a Visual Studio project.

Once the previous command completede successfully, compile the examples and install the **castor** headers:

.. code:: text

    make
    make install

The executable files for the examples can be found in the ``castor/build/demo/demo_*`` subfolders. The folder containing the headers is copied in the ``/mingw64/include/`` sub-directory.

**Remark:** if you have questions or remarks about the installation procedure on Windows, please contact Marc Bakry (contact info at :ref:`label-developpers`).
