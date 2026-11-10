# System-specific makefile for Windows using MSYS2 (UCRT64 environment).
#
# Required packages (run in an MSYS2 UCRT64 shell):
#   pacman -S make git mingw-w64-ucrt-x86_64-{gcc-fortran,msmpi,openblas,scalapack,fftw,libxc}
# The MS-MPI runtime (mpiexec, smpd) must be installed separately from
# Microsoft (e.g. "winget install Microsoft.msmpi").
#
# Build from an MSYS2 UCRT64 shell with:  make SYSTEM=msys2
# Full instructions: "Installing on Windows" in docs/installing.rst.

# Set compilers (MS-MPI wrapper shipped by mingw-w64-ucrt-x86_64-msmpi)
FC=mpifort

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS= -fopenmp
# Set this to "OMP_DUMMY = DUMMY" if compiling without openmp
# Set this to "OMP_DUMMY = " if compiling with openmp
OMP_DUMMY =

# Set BLAS and LAPACK libraries (OpenBLAS provides both)
BLAS= -lopenblas
# Full scalapack library call; remove -lscalapack if using dummy diag module.
SCALAPACK = -lscalapack

# LibXC compatibility (v5, v6 and v7 have the same interface)
XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf03 -lxc
XC_COMPFLAGS = -I/ucrt64/include/libxc

# Set FFT library
FFT_LIB=-lfftw3
FFT_OBJ=fft_fftw3.o

# Set ELPA library
ELPA_LIB =
ELPA_INC =

LIBS= $(FFT_LIB) $(ELPA_LIB) $(XC_LIB) $(SCALAPACK) $(BLAS)

# Compilation flags
# NB for gcc10+ you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 $(OMPFLAGS) $(XC_COMPFLAGS) $(ELPA_INC) -fallow-argument-mismatch

# Linking flags
LINKFLAGS= -L/ucrt64/lib $(OMPFLAGS)

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
# Use dummy ELPAModule or not
ELPA_DUMMY =DUMMY
