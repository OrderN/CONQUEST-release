.. _install:

============
Installation
============

You will need to download and compile the code before you can use it;
we do not supply binaries.

.. _install_down:

Downloading
-----------

CONQUEST is accessed from `the GitHub repository
<https://github.com/OrderN/CONQUEST-release/>`_;
it can be cloned:

``git clone https://github.com/OrderN/CONQUEST-release destination-directory``

where ``destination-directory`` should be set by the user.
Alternatively, it can be downloaded from GitHub as a zip file and
unpacked:

`<https://github.com/OrderN/CONQUEST-release/archive/master.zip>`_

Go to :ref:`top <install>`

.. _install_compile:

Compiling
---------

Once you have the distribution, you will need to compile the main
Conquest code (found in the ``src/`` directory), along with the ion file
generation code (found in the ``tools/`` directory).  Conquest requires
a working MPI installation including a Fortran90 compiler (often
``mpif90`` but this can vary), along with a few standard libraries:

* BLAS and LAPACK (normally provided by the system vendor)
* FFTW 3.x (more detail can be found at `http://www.fftw.org/ <http://www.fftw.org/>`_)
* ScaLAPACK (often provided as part of an HPC system; the source code
  can be obtained from `the netlib repository <http://www.netlib.org/scalapack/>`_ if
  you need to compile it)

Additionally, Conquest can use LibXC if it is available (v4.x or
later).

The library locations are set in the ``system.make`` file in the ``src/system``
directory, along with other parameters needed for compilation.  The default file
name is ``system.make`` but you can select another file with ``make SYSTEM=label``
which would then use the file ``system.label.make`` in the ``src/system`` directory.
``system.<systemname>.make``
files are provided for some HPC systems used by the community, but if you want to run
locally or on a different system, you will need to create an appropriate ``system.make``
file. Use ``src/system/system.example.make`` as a starting point.

* ``FC`` (typically ``FC=mpif90`` will be all that is required)
* ``COMPFLAGS`` (set these to specify compiler options such as
  optimisation)
* ``BLAS`` (specify the BLAS and LAPACK libraries)
* ``SCALAPACK`` (specify the ScaLAPACK library)
* ``FFT_LIB`` (must be left as FFTW)
* ``XC_LIBRARY`` (choose ``XC_LIBRARY=CQ`` for the internal Conquest
  library, otherwise ``XC_LIBRARY=LibXC_v4`` for LibXC v4.x, or ``XC_LIBRARY=LibXC_v5``
  for LibXC v5.x and v6.x)
* Two further options need to be set for LibXC:

  + ``XC_LIB`` (specify the XC libraries)
  + ``XC_COMPFLAGS`` (specify the location of the LibXC include and
    module files, e.g. ``-I/usr/local/include``)

Once these are set, you should make the executable using ``make``.

The ion file generation code is compiled using the same options
required for the main code.

Go to :ref:`top <install>`

Multi-threading
~~~~~~~~~~~~~~~

CONQUEST can use OpenMP for multi-threading; some multi-threading is available throughout the code, while there are specific matrix multiplication routines which can use multi-threading for the linear scaling solver.  The number of threads is set via the environment variable ``OMP_NUM_THREADS``.

Compiler flags to enable OpenMP are dependent on the vendor, but should be specified via ``OMPFLAGS`` in the ``system.make`` file.  If compiling with OpenMP then you should also change the variable ``OMP_DUMMY`` in the same file to be blank to enable the number of threads to be included in the output.

On some systems, the default stack size for OpenMP is set to be rather small, and this can cause a segmentation fault when running with multiple threads.  We recommend testing the effect of the environment variable ``OMP_STACKSIZE`` (and suggest setting it to 50M or larger as a first test).

Most OpenMP multi-threading in CONQUEST uses the ``runtime`` schedule. This means the type of scheduling of work to threads can be set by the user by setting the ``OMP_SCHEDULE`` `variable<https://www.openmp.org/spec-html/5.0/openmpse49.html>`_. If the variable is unset, OpenMP will use a default implementation defined schedule. 

Go to :ref:`top <install>`

.. _install_spack:

Installing with Spack
-----------

CONQUEST and all of its dependencies can be installed with `Spack <https://spack.io/>`_.
The CONQUEST package requires Spack v0.21 or later. If Spack isn't available or up to date on your
system, it is relatively straightforward to install it with user permissions following the
`install instructions <https://spack.readthedocs.io/en/latest/getting_started.html#installation>`_.
When setting up Spack on a new system, it is recommended to configure it to use available
`system compilers <https://spack.readthedocs.io/en/latest/getting_started.html#compiler-configuration>`_
and `system packages <https://spack.readthedocs.io/en/latest/getting_started.html#system-packages>`_.
Once spack is installed and set up, install CONQUEST with:

``spack install conquest``

and load the ``Conquest`` executable to ``PATH`` with

``spack load conquest``

The build can be customized by adding options to the
`Spack spec <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_ ``conquest``.
The CONQUEST package includes variants for OpenMP support and different matrix multiplication kernels; more details can be found in the `Spack CONQUEST package <https://spack.readthedocs.io/en/latest/package_list.html#conquest>`_.

Installing on Ubuntu
-----------

CONQUEST can be compiled on Ubuntu after installing the required packages. The below instructions are given for Ubuntu 22.04 LTS and Ubuntu 24.04 LTS.
The source files will be downloaded into the ``${USER}/local/src`` directory. The ${USER} variable will be automatically replaced by the current username.

Install needed packages
~~~~~~~~~~~~~~~~

.. code-block:: bash

    sudo apt update
    sudo apt upgrade

    sudo apt install -y build-essential                                    # GCC and other tools for software development
    sudo apt install -y openmpi-bin libopenmpi-dev                         # MPI
    sudo apt install -y libfftw3-dev                                       # FFT
    sudo apt install -y libblas-dev liblapack-dev libscalapack-openmpi-dev # Linear algebra

Install libxc
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    cd $HOME && mkdir local
    cd $HOME/local && mkdir src && cd src

    cd $HOME/local/src
    wget https://gitlab.com/libxc/libxc/-/archive/6.2.2/libxc-6.2.2.tar.bz2 -O libxc.tar.bz2
    tar -xf libxc.tar.bz2
    cd libxc-6.2.2 && autoreconf -i && ./configure --prefix=$HOME/local
    make
    make check && make install

Download CONQUEST
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    cd $HOME/local/src
    git clone https://github.com/OrderN/CONQUEST-release.git conquest_master
    cd conquest_master/src

Prepare makefile
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Prepare system.make file for Ubuntu.
    cat > system/system.make << EOF

    # Set compilers
    FC=mpif90
    F77=mpif77

    # Linking flags
    LINKFLAGS= -L\${HOME}/local/lib -L/usr/local/lib -fopenmp
    ARFLAGS=

    # Compilation flags
    # NB for gcc10 you need to add -fallow-argument-mismatch
    COMPFLAGS= -O3 \$(XC_COMPFLAGS) -fallow-argument-mismatch
    COMPFLAGS_F77= \$(COMPFLAGS)

    # Set BLAS and LAPACK libraries
    # Generic
    BLAS= -llapack -lblas

    # Full library call; remove scalapack if using dummy diag module
    LIBS= \$(FFT_LIB) \$(XC_LIB) -lscalapack-openmpi \$(BLAS)

    # LibXC compatibility (LibXC below) or Conquest XC library

    # LibXC compatibility
    # Choose LibXC version: v4 (deprecated) or v5/6 (v5 and v6 have the same interface)
    XC_LIBRARY = LibXC_v5
    XC_LIB = -lxcf90 -lxc
    XC_COMPFLAGS = -I\${HOME}/local/include -I/usr/local/include

    # Set FFT library
    FFT_LIB=-lfftw3
    FFT_OBJ=fft_fftw3.o

    # Matrix multiplication kernel type
    MULT_KERN = default
    # Use dummy DiagModule or not
    DIAG_DUMMY =

    EOF

Compile CONQUEST
~~~~~~~~~~~~~~~

.. code-block:: bash


    dos2unix ./makedeps       # For Windows Subsystem for Linux (WSL), there may be some incompatibilities thus file conversion is recommended.
    make                      # Or make -j`nproc` for parallel compilation using all available cores
Go to :ref:`top <install>`

