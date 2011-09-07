# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: system.make,v 1.2 2010/06/19 22:07:21 lt Exp $
#
# Conquest example system.make
#
#
FC = ftn
F77 = ftn
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS = 
# COMPFLAGS = -g
COMPFLAGS = -fastsse
# COMPFLAGS_F77 = -g
COMPFLAGS_F77 = -fastsse
ARFLAGS = 

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
BLACS = 
LIBS = $(FFT) 
# This line is for systems with NO ScaLAPACK and BLAS - replace DiagModule.f90 with DiagModule.f90_DUMMY
#LIBS = $(FFT) -llapack -lblas $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC=$(F77)" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")
