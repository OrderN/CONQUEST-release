# This is an example system-specific makefile. You will need to adjust
# it for the actual system you are running on.

# Set compilers
FC=mpif90
F77=mpif77

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS= 

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 $(OMPFLAGS) $(XC_COMPFLAGS)
COMPFLAGS_F77= $(COMPFLAGS)

# Set BLAS and LAPACK libraries
# MacOS X
# BLAS= -lvecLibFort
# Intel MKL use the Intel tool
# Generic
BLAS= -llapack -lblas

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
XC_COMPFLAGS = -I/usr/local/include

# Set FFT library
FFT_LIB=-lfftw3
FFT_OBJ=fft_fftw3.o

# Full library call; remove -lscalapack if using dummy diag module.
# If using OpenMPI, use -lscalapack-openmpi instead.
LIBS= $(FFT_LIB) $(XC_LIB) -lscalapack $(BLAS)

# Linking flags
LINKFLAGS= -L/usr/local/lib $(OMPFLAGS)
ARFLAGS=

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
# Use dummy omp_module or not.
# Set this to "OMP_DUMMY = DUMMY" if compiling without openmp
# Set this to "OMP_DUMMY = " if compiling with openmp
OMP_DUMMY = DUMMY

