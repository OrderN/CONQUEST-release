# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: irix_mips4.make,v 1.1 2003/03/12 14:52:22 drb Exp $
#
# Conquest system file for IRIX on mips4
#
#
FC = f90
LINKFLAGS = -64 -mips4
COMPFLAGS = -64 -extend_source -mips4 -O2

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
LIBS = -lcomplib.sgimath -lscalapack -lmpiblacs -lsma -lmpi
# This line is for systems with dummy DiagModule
#LIBS = -lcomplib.sgimath -lsma -lmpi

# Matrix multiplication kernel type
MULT_KERN = default
# Use dummy DiagModule or not
DIAG_DUMMY =

