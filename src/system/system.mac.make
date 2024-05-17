# This is an example system-specific makefile. You will need to adjust
# it for the actual system you are running on.

# Set compilers
FC=/opt/homebrew/bin/mpifort

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS= -fopenmp
# Set this to "OMP_DUMMY = DUMMY" if compiling without openmp
# Set this to "OMP_DUMMY = " if compiling with openmp
OMP_DUMMY = 

# Set BLAS and LAPACK libraries
# MacOS X
# BLAS= -lvecLibFort
# Intel MKL use the Intel tool
# Generic
BLAS= -llapack -lblas
# Full scalapack library call; remove -lscalapack if using dummy diag module.
# If using OpenMPI, use -lscalapack-openmpi instead.
# If using Cray-libsci, use -llibsci_cray_mpi instead.
SCALAPACK = -lscalapack

# LibXC: choose between LibXC compatibility below or Conquest XC library

# Conquest XC library
#XC_LIBRARY = CQ
#XC_LIB =
#XC_COMPFLAGS =

# LibXC compatibility
# Choose LibXC version: v4 (deprecated) or v5/6 (v5 and v6 have the same interface)
#XC_LIBRARY = LibXC_v4
XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf90 -lxc
XC_COMPFLAGS = -I/opt/homebrew/Cellar/libxc/6.2.2/include

# Set FFT library
FFT_LIB=-lfftw3
FFT_OBJ=fft_fftw3.o

LIBS= $(FFT_LIB) $(XC_LIB) $(SCALAPACK) $(BLAS)

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -fallow-argument-mismatch -O3 $(OMPFLAGS) $(XC_COMPFLAGS) -I/opt/homebrew/Cellar/openblas/0.3.27/include -I/opt/homebrew/Cellar/lapack/3.12.0/include -I/opt/homebrew/Cellar/fftw/3.3.10_1/include

# Linking flags
LINKFLAGS= $(OMPFLAGS) -L/opt/homebrew/Cellar/openblas/0.3.27/lib -L/opt/homebrew/Cellar/lapack/3.12.0/lib -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib -L/opt/homebrew/Cellar/libxc/6.2.2/lib -L/opt/homebrew/Cellar/scalapack/2.2.0_1/lib

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
