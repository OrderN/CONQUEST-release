# system.waltz.make for Ubuntu 24.04 + CUDA 12.9 + OpenMPI + ELPA-GPU
FC = mpifort

# Set FFT library
FFT_LIB = -L/usr/lib/x86_64-linux-gnu -lfftw3
FFT_OBJ = fft_fftw3.o

# Set ScaLAPACK and BLAS
SCALAPACK = -L/usr/lib/x86_64-linux-gnu -lscalapack-openmpi
BLAS = -L/usr/lib/x86_64-linux-gnu -llapack -lblas

# Set ELPA library with CUDA GPU support
ELPA_LIB = -L/home/augustin/local/elpa-gpu/lib -lelpa_openmp -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusolver -lstdc++
ELPA_INC = -I/home/augustin/local/elpa-gpu/include/elpa_openmp-2024.05.001/modules
ELPA_DUMMY =

# XC library (native CQ)
XC_LIBRARY = CQ
XC_LIB =
XC_COMPFLAGS =

# OpenMP
OMPFLAGS = -fopenmp

# Flags
COMPFLAGS = -O3 $(OMPFLAGS) -fallow-argument-mismatch $(XC_COMPFLAGS) $(ELPA_INC)
LINKFLAGS = $(OMPFLAGS) -Wl,-rpath,/home/augustin/local/elpa-gpu/lib -Wl,-rpath,/usr/local/cuda/lib64

LIBS = $(FFT_LIB) $(ELPA_LIB) $(XC_LIB) $(SCALAPACK) $(BLAS)
