# This is a system.make file for the UCL Myriad machine. See
# https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/ for details
# module load gcc-libs/10.2.0
# module load compilers/intel/2022.2
# module load mpi/intel/2021.6.0
# module load libxc/6.2.2/intel-2022

# Set compilers
FC=mpif90
F77=mpif77

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS= -fopenmp

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
# Choose LibXC version: v4 (deprecated) or v5 (v5, v6 and v7 have the same interface)
#XC_LIBRARY = LibXC_v4
#XC_LIB = -lxcf90 -lxc
XC_LIBRARY = LibXC_v5
XC_LIB = -lxcf03 -lxc
XC_COMPFLAGS = -I/usr/local/include

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 -g $(OMPFLAGS) $(XC_COMPFLAGS) $(ELPA_INC) -I"${MKLROOT}/include"
COMPFLAGS_F77= $(COMPFLAGS)

# Set FFT library
FFT_LIB=-lmkl_rt
FFT_OBJ=fft_fftw3.o

# Set ELPA library
#ELPA_LIB = -L/**/lib -lelpa
#ELPA_INC = -I/**/modules/
ELPA_LIB =
ELPA_INC =

# Full library call; remove scalapack if using dummy diag module
# If using OpenMPI, use -lscalapack-openmpi instead.
#LIBS= $(FFT_LIB) $(XC_LIB) -lscalapack $(BLAS)
LIBS= $(FFT_LIB) $(ELPA_LIB) $(XC_LIB)

# Linking flags
LINKFLAGS= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -ldl $(OMPFLAGS) $(XC_LIB)
ARFLAGS=

# Matrix multiplication kernel type
MULT_KERN = ompGemm_m
# Use dummy DiagModule or not
DIAG_DUMMY =
# Use dummy omp_module or not.
# Set this to "OMP_DUMMY = DUMMY" if compiling without openmp
# Set this to "OMP_DUMMY = " if compiling with openmp
OMP_DUMMY = 
# Use dummy ELPAModule or not
ELPA_DUMMY =DUMMY

