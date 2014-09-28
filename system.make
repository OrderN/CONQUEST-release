# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: system.make,v 1.2 2013/04/10 22:07:21 lat Exp $
#
# ISM RATM system.make
#
#

MPICHROOT=/opt/mpich2/1.4.1p1/intel11.1
MKLROOT=/opt/intel/Compiler/11.1/072/mkl
FFTWROOT=/opt/fftw/intel11.1

FC=$(MPICHROOT)/bin/mpif90
F77=$(MPICHROOT)/bin/mpif77

# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS= -L/opt/mpich2/1.4.1p1/intel11.1/lib
ARFLAGS=

# L.A. Truflandier intel/mkl/mpich2 compilation
COMPFLAGS= -O1 -fp-model strict -g -traceback -I$(MPICHROOT)/include # -I$(FFTWROOT)/include
COMPFLAGS_F77= $(COMPFLAGS)
LIBS=  $(FFT_LIB)
LIBS+= $(MKLROOT)/lib/em64t/libmkl_scalapack_lp64.a
LIBS+= -Wl,--start-group $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a
LIBS+= $(MKLROOT)/lib/em64t/libmkl_sequential.a
LIBS+= $(MKLROOT)/lib/em64t/libmkl_core.a
LIBS+= $(MKLROOT)/lib/em64t/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
LIBS+= -lpthread

# **<lat>** EXX is "stable" with FFTW (need to fix for the others) 
FFT_LIB=-L$(FFTWROOT)/lib -lfftw3
FFT_OBJ=fft_fftw3.o

#FFT_LIB=-LFFT.FFTE -lffte
#FFT_OBJ=fft_ffte.o
#fft_ffte.o: FFT.FFTE/libffte.a
#FFT.FFTE/libffte.a:
#	(cd FFT.FFTE; $(MAKE) "FC77=$(F77)" "FC90=${FC}" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

#FFT_LIB=-LFFT.GPFA -lgpfa
#FFT_OBJ=fft_gpfa.o
#fft_gpfa.o: FFT.GPFA/libgpfa.a
#FFT.GPFA/libgpfa.a:
#	(cd FFT.GPFA; $(MAKE) "FC77=$(F77)" "FC90=${FC}" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")


# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =


