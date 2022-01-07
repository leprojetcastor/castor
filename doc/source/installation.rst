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

**Note** : on macOS it is recommended to install these dependencies with `brew <https://brew.sh/>`_.

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

Windows (10)
++++++++++++

This solution has been tested on Windows 10 only (but may work on other version). Since there is no *built_in* available package manager, the different components will be installed *by-hands* using only *Windows-like* tools.

C++ compiler and CMake
----------------------

The first step is to install a suitable C++ compiler. In the following instructions we will only use the compiler provided with `Visual Studio <https://visualstudio.microsoft.com/fr/downloads/>`_, version 2019 or later (previous version may work but have not been tested). The Visual Studio framework also provides a customized command prompt named `x64 Native Tools Command Prompt`.

We will also install the `CMake <https://cmake.org>`_ tools. Download the latest binary distribution for Windows. After installing CMake, open the `x64 Native Tools Command Prompt` and execute the command `cmake-gui`. If the command fails, find the install folder of CMake and add the `CMakeInstallFolder\bin` subfolder to the Windows `%PATH%` environment variable, restart the command prompt and try again. The graphical interface of CMake should open (close it for now).

BLAS/LAPACK library
-------------------

The simplest way is probably to install `OpenBLAS <https://www.openblas.net>`_ which implements both interfaces. Compiling the library can quickly become painful as a Fortran compiler is required. Thankfully, the developers have made precompiled binaries available. Installing OpenBLAS can be done following these steps:

1. Go to `https://github.com/xianyi/OpenBLAS/releases <https://github.com/xianyi/OpenBLAS/releases>`_, look for an archive named **OpenBLAS-0.x.x-x64.zip** (or **OpenBLAS-0.x.x-x86.zip** for older architectures) in the section **Assets** corresponding to the version you wish to use and download it. The demos were tested originally with `OpenBLAS 0.3.12` so any later version should be fine.
2. Extract the downloaded archive in a folder of your choice (for example, create a folder `openblas`). This folder should now contain three subdirectories:
    - `openblas\bin` should contain a file named `libopenblas.dll`.
    - `openblas\include` should contain the header files, including `cblas.h`.
    - `openblas\lib` should contain `.lib` files (including `libopenblas.lib`) and a `openblas\lib\cmake` subfolder.
3. Add the subfolder `openblas\bin` to Windows environment variable `%PATH%`.

The BLAS and LAPACK are now ready to use.

**Remark** It is also possible to download the Intel MKL library through the framework https://github.com/oneapi-src/oneMKL or the `Intel website <https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html>`_. However, this implementation features different header names and requires a modification of the source files of **castor** (namely replace `cblas.h` by `mkl_cblas.h` wherever it appears). For this reason, we do not insist further.

VTK framework
-------------

Unfortunately, the developers of the VTK framework do not provide *ready-to-use* binaries meaning that we must compile the sources by ourselves. It is performed as follows:

1. Download the sources of VTK on the main website `https://vtk.org <https://vtk.org>`_. Choose a version of the `9.x.x` branch. Uncompress the archive in a folder of your choice.
2. Open the `x64 Native Tools Command Prompt` and move to the newly created VTK folder (use the `dir pathToFolder` command). Create a *build* folder using `mkdir build` and move to this folder.
3. Execute `cmake-gui ..` which should open the CMake graphical interface. Click on `Configure`, choose the `ninja` generator and keep the default configuration. Finally, click on `Generate`. CMake will generate the build files.
4. Go back to the command prompt and execute the command `ninja`. The compilation of VTK begins and *may* take some time (a few minutes to a few dozen of minutes depending on the computer).
5. Once the compilation is over, execute `ninja install` which will install the library in the default directory `C:\Program Files (x86)\VTK`.
6. The last step is to add the subfolder `C:\Program Files (x86)\VTK\bin` to the Windows `%PATH%` environment variable.

The installation of VTK is now completed.

Compile the demos
-----------------

In this section, we will give the instructions on how to compile the examples of castor. The steps are the following:

1. Download the sources of **castor** from the `main repo <https://gitlab.labos.polytechnique.fr/leprojetcastor/castor.git>`_.
2. Open the `x64 Native Tools Command Prompt` and got to the `castor` folder. Create a `castor\build` directory and move to it.
3. Execute `cmake-gui ..` and click on the `Configure` button. Choose the `ninja` generator on the list and let all other options by default. This last operation *should fail* as CMake cannot find BLAS/LAPACK nor VTK.
4. In the list of CMake variables, look for `VTK_DIR` and set it to the `VTK\lib\cmake\vtk-9.x` folder. Do the same for the BLAS-related variables. Look for the variable `CBLAS_INCLUDE_DIR` and set it to the `openblas\include`subfolder.
5. Click again on `Configure` then on `Generate`.
6. Finally, execute `ninja` in the command prompt to start building the demo executable. The corresponding file can then be found in the `castor\build\demo` subfolder.
