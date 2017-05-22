# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: system.make,v 1.2 2013/04/10 22:07:21 lat Exp $
#
# ISM RATM system.make
#
#

FC=mpif90
F77=mpif77

# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS= -L/opt/local/lib
ARFLAGS=

# L.A. Truflandier intel/mkl/mpich2 compilation
COMPFLAGS= -g -O3 -fbounds-check
COMPFLAGS_F77= $(COMPFLAGS)
LIBS=  $(FFT_LIB) -lscalapack ${BLACS} -latlas -llapack -lf77blas -lcblas

# **<lat>** EXX is "stable" with FFTW (need to fix for the others) 
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


