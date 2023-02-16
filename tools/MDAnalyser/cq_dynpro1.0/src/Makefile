# Makefile for CP/Dynpro from QEspresso v5.1d
# Mon May 30 19:00:53 CEST 2016

include	../make.sys

# location of needed modules and included files (if any)
MODFLAGS = $(MOD_FLAG)../modules 

IFLAGS=

DYNMODS = \
./dynpro_math.o \
./dynpro_inp.o \
./dynpro_cell.o \
./dynpro_set.o \
./dynpro_ran.o \
./dynpro_cqc.o \
./dynpro_rdf.o \
./dynpro_fft.o \
./dynpro_psd.o \
./dynpro_mod.o \
./dynpro_sub.o \
./dynpro_sph.o 

QEMODS = ../modules/libqemod.a

LIBOBJS = ../flib/flib.a 


all : dynpro_version dynpro.x

dynpro_version : dynpro_ver.h

dynpro_ver.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> dynpro_ver.h

dynpro.x : $(DYNMODS) dynpro.o $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o dynpro.x dynpro.o $(DYNMODS) $(QEMODS) $(LIBOBJS) 
	- (cd ../bin ; ln -fs ../src/dynpro.x . )


libdynpro.a : $(FOBJS) $(LOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f dynpro.x *.o *.mod dynprov.h *.i core* *.F90 fort* \
		*.cpp *.d *.L *.a *.s 
	-/bin/rm -rf ../bin/dynpro.x

#include make.depend
