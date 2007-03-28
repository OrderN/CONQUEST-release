# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: irix_mips4.make,v 1.1 2003/03/12 14:52:22 drb Exp $
#
# Conquest system file for IRIX on mips4
#
#
FC = mpif77
F77 = mpif77
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS = 
COMPFLAGS = -O2 -check all -fltconsistency -fpe0 -traceback
COMPFLAGS_F77 = -O2 -check all -fltconsistency -fpe0 -traceback
ARFLAGS = 

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
LIBS = $(FFT) -lscalapack -lblacsF77init -lblacs -L/opt/intel/mkl721/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lmpi $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC=$(F77)" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")

#Input library
FDF=libfdf.a
$(FDF):
	(cd fdf; $(MAKE) "FC=$(FC)" "F77=$(F77)" "FFLAGS=-g " "ARFL=$(ARFLAGS)" module)
