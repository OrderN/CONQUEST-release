#

# Set compilers
FC=mpif90
F77=mpif77

# OpenMP flags
OMPFLAGS= -fopenmp

# Linking flags
LINKFLAGS=  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -ldl $(OMPFLAGS)
ARFLAGS=

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 $(OMPFLAGS) $(XC_COMPFLAGS) -I"${MKLROOT}/include"
COMPFLAGS_F77= $(COMPFLAGS)

# Full library call; remove scalapack if using dummy diag module
LIBS= $(FFT_LIB) $(XC_LIB)

# LibXC compatibility (LibXC below) or Conquest XC library

# LibXC compatibility
# Choose LibXC version: v4 or v5
XC_LIBRARY = LibXC_v4
# XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf90 -lxc
XC_COMPFLAGS = -I/shared/ucl/apps/libxc/4.2.3/intel-2018/include

# Set FFT library
FFT_LIB=-lmkl_rt
FFT_OBJ=fft_fftw3.o

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =


