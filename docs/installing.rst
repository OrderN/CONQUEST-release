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

Additionally, Conquest can use LibXC if it is available (v2.x or
later).

The library locations are set in the ``system.<systemname>.make`` file in the ``src/system``
directory, along with other parameters needed for compilation. ``system.<systemname>.make``
files are provided for some HPC systems used by the community, but if you want to run
locally or on a different system, you need to provide an appropriate ``system.<systemname>.make``
file. Use ``src/system/system.example.make`` as a starting point. Get the ``<systemname>``
by running ``hostname -d`` in your prompt, then name your file appropriately and move it to
the ``src/system`` directory. If ``hostname -d`` returns empty (e.g. you are running on a
local machine), the system-specific makefile should be named ``system.make``.

* ``FC`` (typically ``FC=mpif90`` will be all that is required)
* ``COMPFLAGS`` (set these to specify compiler options such as
  optimisation)
* ``BLAS`` (specify the BLAS and LAPACK libraries)
* ``SCALAPACK`` (specify the ScaLAPACK library)
* ``FFT_LIB`` (must be left as FFTW)
* ``XC_LIBRARY`` (choose ``XC_LIBRARY=CQ`` for the internal Conquest
  library, otherwise ``XC_LIBRARY=LibXC_v2or3`` for LibXC v2.x or v3.x, or ``XC_LIBRARY=LibXC_v4``
  for LibXC v4.x)
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

Compiler flags to enable OpenMP are dependent on the vendor, but should be specified via ``COMPFLAGS`` and ``LINKFLAGS`` in the ``system.make`` file.  If compiling with OpenMP then you should also change the variable ``OMP_DUMMY`` in the same file to be blank to enable the number of threads to be included in the output.

On some systems, the default stack size for OpenMP is set to be rather small, and this can cause a segmentation fault when running with multiple threads.  We recommend testing the effect of the environment variable ``OMP_STACKSIZE`` (and suggest setting it to 50M or larger as a first test).

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

Go to :ref:`top <install>`

