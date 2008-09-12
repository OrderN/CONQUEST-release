# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id$
#
# Conquest example system.make
#
#
FC = mpif90
F77 = mpif77
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS = 
COMPFLAGS = -O2 
COMPFLAGS_F77 = -O2 
ARFLAGS = 

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
LIBS = $(FFT) -L/scratch/drb/lib -lscalapack -lmpiblacsF77init -lmpiblacs -lacml $(FFT)
# This line is for systems with NO ScaLAPACK and BLAS - replace DiagModule.f90 with DiagModule.f90_DUMMY
#LIBS = $(FFT) -llapack -lblas $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC=$(F77)" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")
