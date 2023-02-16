# Makefile for flib

include ../make.sys

OBJS = \
inpfile.o \
int_to_char.o \

all: flib.a 

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@


clean :
	- /bin/rm -f *.a *.o *.mod *.i *.F90 core* *.L

#include make.depend
