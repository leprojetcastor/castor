.. _label-installation:

Installation
============

The simplest way to get the dense matrix part of the framework **castor** is to download the `latest version of header files <https://gitlab.labos.polytechnique.fr/leprojetcastor/castor/-/jobs/artifacts/master/download?job=deploy>`_ and include `matrix.hpp` during the compilation of your c++ program using the library, see :ref:`label-compilation` in the :ref:`label-basic` section.

For a complete installation integrating check dependencies and examples compilation, we describe below the procedure.

**MacOS** (11 and 12) and **Ubuntu** (20.04)
++++++++++++++++++++++++++++++++++++++++++++

Installing the dependencies
---------------------------

The **castor** framework depends on two external dependencies : a BLAS/LAPACK implementation in order to use optimized linear algebra, and VTK for the graphical rendering.

The BLAS/LAPACK implementation which has been tested are :

- MKL 2020.0.166 : `MKL informations <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>`_,   
- OpenBlas 0.3.19 : `OpenBLAS informations <https://www.openblas.net/>`_,   
- vecLib : from accelerate framework in MacOS.   

The version of VTK library which has been tested is `9.1.0 <https://vtk.org/download/>`_.

The last tool to perform installation is ``CMake``, at least version `3.18  <https://cmake.org/download/>`_.

Note : on macOS it is recommended to install these dependencies with `brew <https://brew.sh/>_`

From git repository with CMake
------------------------------

You can install the **castor** library from source with ``CMake``. The source files of the library is available here `<https://gitlab.labos.polytechnique.fr/leprojetcastor/castor>`_.

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

Xcode project for macOS
-----------------------

To create a Xcode project with ``CMake``, add option `-G Xcode` in the `cmake` command :
 
.. code::

    $ cmake -G Xcode -DCMAKE_INSTALL_PREFIX=path/to/install/directory 

This project is created in the `build` directory.
