# Compile flags
FC?=mpif77
F90?=mpif90
CFLAGS=
FFLAGS?= -fast $(CFLAGS)
LDFLAGS?= $(CFLAGS)
ARFL?=

OBJS =	gpfa.o setgpfa.o gpfa2f.o gpfa3f.o gpfa5f.o gpfft3.o

libgpfa.a: $(OBJS)
	ar $(ARFL) rv libgpfa.a $(OBJS)
	cp libgpfa.a ..

.f.o:	
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(F90) $(FFLAGS) -c $<

clean:
	rm -f *.o *.a
