# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: system.make,v 1.2 2010/06/19 22:07:21 lt Exp $
#
# LCN Steinway system.make
#
#
FC = mpif90
F77 = mpif77
# If you need 64-bit flags for compiling, add them for linking
LINKFLAGS = 
COMPFLAGS = -fast -tp k8-64
# COMPFLAGS = -g
COMPFLAGS_F77 = -fast -tp k8-64
# COMPFLAGS_F77 = -g
ARFLAGS = 

# This line is for systems with ScaLAPACK, BLACS and diagonalisation
BLACS = /usr/local/pgi/linux86-64/2010/mpi/mpich/lib/blacsCinit_MPI-LINUX-0.a /usr/local/pgi/linux86-64/2010/mpi/mpich/lib/blacs_MPI-LINUX-0.a /usr/local/pgi/linux86-64/2010/mpi/mpich/lib/blacsF77init_MPI-LINUX-0.a 
LIBS = -L/usr/local/pgi/linux86-64/10.0/libso/ $(FFT) -lscalapack ${BLACS} -lacml -lacml_mv $(FFT) 
# This line is for systems with NO ScaLAPACK and BLAS - replace DiagModule.f90 with DiagModule.f90_DUMMY
#LIBS = $(FFT) -llapack -lblas $(FFT)

# Default FFT (GPFA) - replace as necessary
FFT=libgpfa.a
$(FFT):
	(cd FFT; $(MAKE) "FC77=$(F77)" "FC90=$(FC)" "FFLAGS=$(COMPFLAGS_F77)" "ARFL=$(ARFLAGS)")
