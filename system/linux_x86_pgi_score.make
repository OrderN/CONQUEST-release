# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: linux_x86_pgi_score.make,v 1.1 2005/05/11 00:41:08 drb Exp $
#
# Conquest system file for PGI on x86 linux with score
#           T. Miyazaki
#
FC = mpif90 -compiler pgi -nochoicemod
LD = mpif90 -compiler pgi -nochoicemod
CC = mpicc
# D.R.Bowler, 2003/07/22
# Added -Bstatic to solve library problems under RedHat 7.3/9
LINKFLAGS = -g77libs  -Bstatic
COMPFLAGS = -g
#COMPFLAGS = -fast
#COMPFLAGS = -fast

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
LIBS = -lscalapack -lblacs -lblacsF77init -lblacs -llapack -lblas $(FFT) dummy.o
# This line is for systems with dummy DiagModule
#LIBS = -L/usr/local/intel/mkl/lib/32/ -lmkl_lapack -lmkl_def -lguide -lpthread $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC=$(FC)" "FFLAGS=$(COMPFLAGS)")

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =