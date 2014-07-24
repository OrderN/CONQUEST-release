# Compile flags
FC=mpif90
CFLAGS=
#FFLAGS= -r8 -fast $(CFLAGS)
LDFLAGS= $(CFLAGS)
ARFL=

OBJS =  gpfa.o setgpfa.o gpfa2f.o gpfa3f.o gpfa5f.o gpfft3.o

libgpfa.a: $(OBJS)
	ar $(ARFL) rv libgpfa.a $(OBJS)
###	cp libgpfa.a ..

.f.o:
	$(FC) $(FFLAGS) -c $<
#.f90.o:
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
#	$(FC) $(COMPFLAGS) -c $<


clean:
	rm -f *.o *.a
