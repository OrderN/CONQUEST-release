# Set compilers
FC=mpifort
F77=mpif77
# OMP flag for compiler and linker
OMPFLAG= -fopenmp
# Linking flags
LINKFLAGS= -L/usr/lib -L/usr/lib/x86_64-linux-gnu $(OMPFLAG)
ARFLAGS=
# Set BLAS and LAPACK libraries
BLAS= -llapack -lblas
# LibXC compatibility (LibXC below) or Conquest XC library
# Conquest XC library
XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf90 -lxc
XC_COMPFLAGS = -I/usr/include
# Set FFT library
FFT_LIB=-lfftw3
FFT_OBJ=fft_fftw3.o
# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
# Full library call; remove scalapack if using dummy diag module
LIBS= $(XC_LIB) -lscalapack-openmpi $(BLAS) $(FFT_LIB)
# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
FPPFLAGS= -DDEBUG
COMPFLAGS= -cpp -g -O3 $(OMPFLAG) $(XC_COMPFLAGS) -fallow-argument-mismatch $(FPPFLAGS)
COMPFLAGS_F77= $(COMPFLAGS)
