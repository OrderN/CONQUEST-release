# This is a system-specific makefile for the AMD aac6 machine used during the
# hackathon in 11-13 March 2026
#
# First load the following modules
# $ module load rocm/7.2.0 fftw/3.3.10 openmpi/5.0.10-ucc1.6.0-ucx1.19.1-xpmem-2.7.4 petsc/3.24.1 elpa/master-c7234ec

# export OMPI_FC=amdflang
# Set compilers
FC=mpif90

# OpenMP flags
# Set this to "OMPFLAGS= " if compiling without openmp
# Set this to "OMPFLAGS= -fopenmp" if compiling with openmp
OMPFLAGS=-fopenmp --offload-arch=gfx942
# OMPFLAGS=-fopenmp 
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
#SCALAPACK =
# SCALAPACK = -lscalapack
# SCALAPACK=-L/shared/apps/ubuntu/rocmplus-afar-22.2.0/petsc/petsc/lib/ -lscalapack-openmpi
SCALAPACK=-L/shared/apps/ubuntu/opt/rocmplus-7.2.0/petsc-v3.24.1/petsc/lib -lscalapack

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
FFT_LIB=-lfftw3
FFT_LIB=-L/shared/apps/ubuntu/opt/rocmplus-7.2.0/fftw-v3.3.10/lib -lfftw3
# FFT_LIB=-L/shared/apps/ubuntu/rocmplus-afar-22.2.0/fftw-v3.3.10/lib -lfftw3
FFT_OBJ=fft_fftw3.o fft_hipfft.o

# Set ELPA library
#ELPA_LIB = -L/shared/apps/ubuntu/opt/rocmplus-7.2.0/elpa-vmaster-c7234ec/lib -lelpa
#ELPA_INC = -I/shared/apps/ubuntu/opt/rocmplus-7.2.0/elpa-vmaster-c7234ec/include/elpa-2026.02.001/modules/
ELPA_LIB = 
ELPA_INC = 
#ELPA_LIB = -L/**/lib -lelpa
#ELPA_INC = -I/**/modules/

LIBS= $(FFT_LIB) $(ELPA_LIB) $(XC_LIB) $(SCALAPACK) $(BLAS)

# Compilation flags
# NB for gcc10 you need to add -fallow-argument-mismatch
COMPFLAGS= -O3 -g $(OMPFLAGS) $(XC_COMPFLAGS) $(ELPA_INC) -I${ROCM_PATH}/include/hipfort/amdgcn

# Linking flags
# LINKFLAGS= -L/usr/local/lib $(OMPFLAGS) -lflang_rt.hostdevice #runtime
LINKFLAGS= -L/usr/local/lib $(OMPFLAGS)  -L${ROCM_PATH}/lib -lhipfort-amdgcn  -lamdhip64 -lhipfft

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
# Use dummy ELPAModule or not
ELPA_DUMMY =DUMMY
#ELPA_DUMMY =
