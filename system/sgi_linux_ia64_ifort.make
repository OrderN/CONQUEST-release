# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: sgi_linux_ia64_ifort.make,v 1.1 2005/07/11 10:11:11 drb Exp $
#
# Conquest system file for IRIX on mips4
#
#
FC = ifort
LINKFLAGS = 
COMPFLAGS = -O2

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
LIBS = -lsdsm -lscs -lmpi $(FFT)
# This line is for systems with dummy DiagModule
#LIBS = -lcomplib.sgimath -lsma -lmpi

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC=$(FC)" "FFLAGS=$(COMPFLAGS)")
