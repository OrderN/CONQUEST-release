# Makefile for Modules

include ../make.sys

# location of needed modules

MODULES = \
kinds.o \
constants.o \
parameters.o \
io_global.o  \
io_files.o \
parallel_include.o \
mp.o \
mp_world.o \
mp_images.o \
mp_pools.o \
mp_pots.o \
mp_bands.o \
mp_diag.o \
command_line_options.o \
mp_global.o 

all : libqemod.a

#all : version version.o libqemod.a

## If no update_version is found, or it is not executable,
## just copy version.f90.in (that contains 'version=unknown').
## We use cat and not cp so that version.f90 has a new timestamp,
## forcing recompilation of version.f90.
## We use the inner if to avoid recompilation at each make.
#version:
#	- ( if test -x ../install/update_version ; then \
#	../install/update_version; \
#	else if test ! -f version.f90 ; then \
#	cat version.f90.in > version.f90 ; fi ; fi )	
#
#version.f90: version

libqemod.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L version.f90

# In this way, the script for generating version.f90 always runs
# (possibly updating the version if you execute "svn update").
# the update_version script should take care of not changing the file 
# if the svn version did not change
.PHONY: version all clean

#include make.depend
