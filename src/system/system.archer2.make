# This is a system-specific makefile for ARCHER2. 
# See https://docs.archer2.ac.uk/ for user documentation
# This works with the gnu programming environment
# Starting from default modules, use
#   module rm PrgEnv-cray
#   module load PrgEnv-gnu

# Set compilers
FC=mpif90

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS=-fopenmp
# Set this to "OMP_DUMMY = " if compiling without openmp
# Set this to "OMP_DUMMY = " if compiling with openmp
OMP_DUMMY= 

# Set BLAS and LAPACK libraries
BLAS=
# Full scalapack library call; remove -lscalapack if using dummy diag module.
# If using OpenMPI, use -lscalapack-openmpi instead.
# If using Cray-libsci, use -llibsci_cray_mpi instead.
SCALAPACK=

# LibXC: choose between LibXC compatibility below or Conquest XC library

# Conquest XC library
XC_LIBRARY = CQ
XC_LIB =
XC_COMPFLAGS =

# LibXC compatibility. Archer2 doesn't provide a module for LibXC, to use it
# install it under your user directory and link as in the example below
#XC_LIBRARY = LibXC_v5
#XC_LIB = -L/mnt/lustre/a2fs-work2/work/ecseah10/ecseah10/tk-ecse08/excalibur-tests/benchmarks/spack/archer2/compute-node/opt/linux-sles15-zen2/gcc-11.2.0/libxc-5.2.3-sq6g4ckyvauupz6rfd2f5uwqp47cb5bl/lib -lxcf90 -lxc
#XC_COMPFLAGS = -I/mnt/lustre/a2fs-work2/work/ecseah10/ecseah10/tk-ecse08/excalibur-tests/benchmarks/spack/archer2/compute-node/opt/linux-sles15-zen2/gcc-11.2.0/libxc-5.2.3-sq6g4ckyvauupz6rfd2f5uwqp47cb5bl/include

# Set FFT library
FFT_LIB=-L/opt/cray/pe/fftw/3.3.10.3/x86_rome/lib -lfftw3
FFT_OBJ=fft_fftw3.o

LIBS= $(FFT_LIB) $(XC_LIB) $(BLAS)

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 -fallow-argument-mismatch -fopenmp  $(XC_COMPFLAGS)

# Linking flags
LINKFLAGS=  -fopenmp -L$(LIBSCI_BASE_DIR)/gnu/9.1/x86_64/lib -lsci_gnu_mpi_mp -lsci_gnu_mp


# Matrix multiplication kernel type
MULT_KERN = ompGemm_m
# Use dummy DiagModule or not
DIAG_DUMMY =
