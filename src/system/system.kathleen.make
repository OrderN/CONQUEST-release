# This is a system.make file for the UCL Kathleen machine. See
# https://www.rc.ucl.ac.uk/docs/Clusters/Kathleen/ for details


# Set compilers
FC=mpif90

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
#BLAS= -llapack -lblas

# LibXC: choose between LibXC compatibility below or Conquest XC library

# Conquest XC library
#XC_LIBRARY = CQ
#XC_LIB =
#XC_COMPFLAGS =

# LibXC compatibility
# Choose LibXC version: v4 (deprecated) or v5/6 (v5 and v6 have the same interface)
# XC_LIBRARY = LibXC_v4
XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf90 -lxc
XC_COMPFLAGS = -I/usr/local/include

# Set FFT library
FFT_LIB=-lmkl_rt
FFT_OBJ=fft_fftw3.o

# Full library call; remove scalapack if using dummy diag module
# If using OpenMPI, use -lscalapack-openmpi instead.
#LIBS= $(FFT_LIB) $(XC_LIB) -lscalapack $(BLAS)
LIBS= $(FFT_LIB) $(XC_LIB)

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -xAVX -O3 -g $(OMPFLAGS) $(XC_COMPFLAGS) -I"${MKLROOT}/include"

# Linking flags
LINKFLAGS= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -ldl $(OMPFLAGS) $(XC_LIB)

# Matrix multiplication kernel type
MULT_KERN = ompGemm
# Use dummy DiagModule or not
DIAG_DUMMY =
