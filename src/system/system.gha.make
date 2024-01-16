# System-specific makefile for the GitHub Actions runners.

# Linking flags
LINKFLAGS=-fopenmp -L/usr/lib -L/usr/lib/x86_64-linux-gnu
ARFLAGS=
# Set BLAS and LAPACK libraries
BLAS= -llapack -lblas
# LibXC compatibility (LibXC below) or Conquest XC library
XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf90 -lxc
XC_COMPFLAGS = -I/usr/include
# Set FFT library
FFT_LIB=-lfftw3
FFT_OBJ=fft_fftw3.o
# Use dummy DiagModule or not
DIAG_DUMMY =
# Full library call; remove scalapack if using dummy diag module
LIBS= $(XC_LIB) -lscalapack-openmpi $(BLAS) $(FFT_LIB)
# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 -fopenmp $(XC_COMPFLAGS) -fallow-argument-mismatch
COMPFLAGS_F77= $(COMPFLAGS)
