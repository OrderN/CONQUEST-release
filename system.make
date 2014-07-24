# -*- mode: makefile -*-
#
# $Id: system.make,v 1.2 2010/06/19 22:07:21 lt Exp $
#
# LCN Steinway system.make
#
#
FC = mpif90
F77 = mpif77
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS = -L${SCALAPACK_LIB}
COMPFLAGS = -O3
# COMPFLAGS = -g
COMPFLAGS_F77 = -O3
# COMPFLAGS_F77 = -g
ARFLAGS =

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
BLACS =
LIBS = $(FFT_LIB) -lscalapack ${BLACS} -llapack -lblas
# This line is for systems with NO ScaLAPACK and BLAS - replace DiagModule.f90 with DiagModule.f90_DUMMY
#LIBS = $(FFT_LIB) -llapack -lblas

#----- Choice of FFT library -----#
### use two lines below if you want to use FFTW3 as FFT library.
### N.B. you need the FORTRAN interface
#FFT_LIB=-L_PUT_LOCATION_OF_FFTW_LIBRARY_HERE_ -lfftw3
#FFT_OBJ=fft_fftw3.o

### use five lines below if you want to use FFTE as FFT library.
#FFT_LIB=-LFFT.FFTE -lffte
#FFT_OBJ=fft_ffte.o
#fft_ffte.o: FFT.FFTE/libffte.a
#FFT.FFTE/libffte.a:
#	(cd FFT.FFTE; $(MAKE) "FC77=$(F77)" "FC90=${FC}" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

### use five lines below if you want to use GPFA as FFT library.
FFT_LIB=-LFFT.GPFA -lgpfa
FFT_OBJ=fft_gpfa.o
fft_gpfa.o: FFT.GPFA/libgpfa.a
FFT.GPFA/libgpfa.a:
	(cd FFT.GPFA; $(MAKE) "FC77=$(F77)" "FC90=${FC}" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =
