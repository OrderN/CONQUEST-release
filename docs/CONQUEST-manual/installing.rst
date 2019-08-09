============
Installation
============

You will need to download and compile the code before you can use it;
we do not supply binaries.

Downloading
-----------

CONQUEST is accessed from `the GitHub repository
<https://github.com/OrderN/CONQUEST-release/>`_;
it can be cloned:

``git clone https://github.com/OrderN/CONQUEST-release destination-directory``

Alternatively, it can be downloaded from GitHub as a zip file and
unpacked:

`<https://github.com/OrderN/CONQUEST-release/archive/master.zip>`_

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

The library locations are set in the ``system.make`` file in the ``src/``
directory, along with other parameters needed for compilation.

* ``FC`` (typically ``FC=mpif90`` will be all that is required)
* ``COMPFLAGS`` (set these to specify compiler options such as
  optimisation)
* ``BLAS`` (specify the BLAS and LAPACK libraries)
* ``SCALAPACK`` (specify the ScaLAPACK library)
* ``FFT_LIB`` (**we may remove this and require FFTW**)
* ``XC_LIBRARY`` (choose ``XC_LIBRARY=CQ`` for the internal Conquest
  library, otherwise ``XC_LIBRARY=LibXC_v2`` or ``XC_LIBRARY=LibXC``
  for LibXC v2.x or higher)
* Two further options need to be set for LibXC:

  + ``XC_LIB`` (specify the XC libraries)
  + ``XC_COMPFLAGS`` (specify the location of the LibXC include and
    module files, e.g. ``-I/usr/local/include``)

Once these are set, you should make the executable using ``make``.

The ion file generation code is compiled using the same options
required for the main code.

Where next?
-----------

While the tutorials have covered the basic operations of Conquest,
there are many more subtle questions and issues.  Detailed
descriptions of the files required and produced by Conquest are given
in :ref:`input-output`, while the input tags themselves are documented
in :ref:`input_tags`.  A set of :ref:`important` give more details
about the operations of the code and the background theory.
