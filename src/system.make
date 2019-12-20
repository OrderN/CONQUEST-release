#

# Set compilers
FC=mpif90
F77=mpif77

# Linking flags
LINKFLAGS= -L/opt/local/lib
ARFLAGS=

# Compilation flags
COMPFLAGS= -O3 $(XC_COMPFLAGS)
COMPFLAGS_F77= $(COMPFLAGS)

# Set BLAS and LAPACK libraries
BLAS= -latlas -llapack -lf77blas -lcblas

# Full library call; remove scalapack if using dummy diag module
LIBS= $(FFT_LIB) $(XC_LIB) -lscalapack $(BLAS)

# LibXC compatibility (LibXC below) or Conquest XC library

# Conquest XC library
#XC_LIBRARY = CQ
#XC_LIB =
#XC_COMPFLAGS =

# LibXC compatibility
# Choose old LibXC (v2.x) or modern versions
#XC_LIBRARY = LibXC_v2
XC_LIBRARY = LibXC
XC_LIB = -lxcf90 -lxc
XC_COMPFLAGS = -I/opt/local/include

# Set FFT library

# EXX is "stable" with FFTW (need to fix for the others) 
FFT_LIB=-lfftw3
FFT_OBJ=fft_fftw3.o

#FFT_LIB=-LFFT.FFTE -lffte
#FFT_OBJ=fft_ffte.o
#fft_ffte.o: FFT.FFTE/libffte.a
#FFT.FFTE/libffte.a:
#	(cd FFT.FFTE; $(MAKE) "FC77=$(F77)" "FC90=${FC}" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

#FFT_LIB=-LFFT.GPFA -lgpfa -L$(FFTWROOT)/lib -lfftw3
#FFT_OBJ=fft_gpfa.o
#fft_gpfa.o: FFT.GPFA/libgpfa.a
#FFT.GPFA/libgpfa.a:
#	(cd FFT.GPFA; $(MAKE) "FC77=$(F77)" "FC90=${FC}" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")


# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =


