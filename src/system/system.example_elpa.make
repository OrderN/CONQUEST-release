# This is an example system-specific makefile. You will need to adjust
# it for the actual system you are running on.

# Set compilers
FC=mpif90

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS=
# Set this to "OMP_DUMMY = DUMMY" if compiling without openmp
# Set this to "OMP_DUMMY = " if compiling with openmp
OMP_DUMMY = DUMMY

# Set BLAS and LAPACK libraries
# MacOS X
# BLAS= -lvecLibFort
# Intel MKL use the Intel tool
# Generic
BLAS= -llapack -lblas
# Full scalapack library call; remove -lscalapack if using dummy diag module.
# If using OpenMPI, use -lscalapack-openmpi instead.
# If using Cray-libsci, use -llibsci_cray_mpi instead.
#SCALAPACK =
#SCALAPACK = -lscalapack
SCALAPACK=-L/shared/apps/ubuntu/rocmplus-7.2.0/petsc/petsc/lib/ -lscalapack-openmpi
#SCALAPACK=-L/shared/apps/ubuntu/rocmplus-7.2.0/petsc/petsc/lib/ -lscalapack

# LibXC: choose between LibXC compatibility below or Conquest XC library

# Conquest XC library
XC_LIBRARY = CQ
XC_LIB =
XC_COMPFLAGS =

# LibXC compatibility
# Choose LibXC version: v4 (deprecated) or v5/6/7 (v5, v6 and v7 have the same interface)
#XC_LIBRARY = LibXC_v4
#XC_LIBRARY = LibXC_v5
#XC_LIB = -lxcf90 -lxc
#XC_LIB = -lxcf03 -lxc
#XC_COMPFLAGS = -I/shared/apps/ubuntu/libxc/include/

# Set FFT library
#FFT_LIB=-lfftw3
FFT_LIB=-L/shared/apps/ubuntu/rocmplus-7.2.0/fftw-v3.3.10/lib -lfftw3
FFT_OBJ=fft_fftw3.o

# Set ELPA library
#ELPA_LIB = -L/shared/apps/ubuntu/rocmplus-7.2.0/elpa/lib -lelpa
#ELPA_INC = -I/shared/apps/ubuntu/rocmplus-7.2.0/elpa/include/elpa-2025.06.001/modules/
ELPA_LIB = 
ELPA_INC = 
#ELPA_LIB = -L/**/lib -lelpa
#ELPA_INC = -I/**/modules/

LIBS= $(FFT_LIB) $(ELPA_LIB) $(XC_LIB) $(SCALAPACK) $(BLAS)

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 $(OMPFLAGS) $(XC_COMPFLAGS) $(ELPA_INC)

# Linking flags
LINKFLAGS= -L/usr/local/lib $(OMPFLAGS)

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
# Use dummy ELPAModule or not
#ELPA_DUMMY =DUMMY
ELPA_DUMMY =
