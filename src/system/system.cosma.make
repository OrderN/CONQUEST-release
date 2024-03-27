# system.make for cosma8 (apr 2023)
# For user docs see https://cosma.readthedocs.io/en/latest/
# load these modules
# module load intel_comp/2022.3.0 compiler mpi mkl
# module load fftw/3.3.10cosma8


# Set compilers
#FC=scorep --user mpif90
#F77=scorep --user mpif77
FC=mpif90
F77=mpif77

# Linking flags
LINKFLAGS=  ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -fopenmp
ARFLAGS=

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -g -fopenmp -O3 $(XC_COMPFLAGS)  -I${MKLROOT}/include/intel64/lp64 -I"${MKLROOT}/include" -fno-omit-frame-pointer -xHost
COMPFLAGS_F77= $(COMPFLAGS)

# LibXC compatibility (LibXC below) or Conquest XC library

# Conquest XC library
XC_LIBRARY = CQ
XC_LIB =
XC_COMPFLAGS =

# Set FFT library
FFT_LIB=-lmkl_rt
FFT_OBJ=fft_fftw3.o

# Full library call; remove scalapack if using dummy diag module
LIBS= $(FFT_LIB) $(XC_LIB)

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =


