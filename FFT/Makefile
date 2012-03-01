# Compile flags
FC77?=mpif77
FC90?=mpif90
CFLAGS=
FFLAGS?= -fast $(CFLAGS)
LDFLAGS?= $(CFLAGS)
ARFL?=

OBJS =	gpfa.o setgpfa.o gpfa2f.o gpfa3f.o gpfa5f.o gpfft3.o

libgpfa.a: $(OBJS)
	ar $(ARFL) rv libgpfa.a $(OBJS)
	cp libgpfa.a ..

.f.o:	
	$(FC77) $(FFLAGS) -c $<
%.o: %.f90
	$(FC90) $(FFLAGS) -c $<

clean:
	rm -f *.o *.a
