# -*- mode:makefile -*-
# $Id: Makefile,v 1.4 2003/03/12 15:04:28 drb Exp $
#
# Makefile for fdf
#
# Adapted for Conquest by D.R.Bowler
FC?=pgf90
F77?=pgf77
CFLAGS= 
FFLAGS?= -O0     $(CFLAGS)
LDFLAGS= $(CFLAGS)
ARFL?=

module:  libfdf.a
libfdf.a:	fdf.o fdf_mod.o parse.o io.o
	ar $(ARFL) qv libfdf.a fdf.o fdf_mod.o parse.o io.o
	cp libfdf.a ..
parse.o:  parse.f 
	$(FC) $(FFLAGS) -c parse.f
fdf_mod.o: fdf_mod.f parse.o
	$(FC) $(FFLAGS) -c fdf_mod.f
# This line below is for "normal" compilers that produce a .mod file
#	@cp -p *.mod ..
# These lines are for the Intel ifc which doesn't !
	@cp -p fdf_mod.f parse.f ..
	(cd ..;	$(FC) $(FFLAGS) -c parse.f ; $(FC) $(FFLAGS) -c fdf_mod.f )
fdf.o:  fdf.f fdf.h
	$(FC) $(FFLAGS) -c fdf.f
#
clean:
	rm -f *.o libfdf.a *.mod fdf_mod.o *.d work.pc work.pcl ../fdf_mod.d ../fdf_mod.o ../parse.d ../parse.o ../libfdf.a
#
.f.o:
	$(FC) $(FFLAGS) -c $<


