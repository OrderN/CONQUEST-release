# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: unicos_cray.make,v 1.1 2003/03/12 14:52:22 drb Exp $
#
# Conquest system file for IRIX on mips4
#
#
FC = f90
LINKFLAGS =
COMPFLAGS = -N132 -Oscalar3,pipeline2,unroll1,bl

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
# Not known how to link
#LIBS = -lblas -lmpi
# This line is for systems with dummy DiagModule
LIBS = -lblas -lmpi

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =