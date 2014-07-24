# Compile flags
FC=mpif90
CFLAGS=
#FFLAGS= -r8 -fast $(CFLAGS)
LDFLAGS= $(CFLAGS)
ARFL=

OBJS =  zfft1d.o zfft3d.o fft235.o kernel.o factor.o

libgpfa.a: $(OBJS)
	ar $(ARFL) rv libffte.a $(OBJS)
###	cp libgpfa.a ..

.f.o:
	$(FC) $(FFLAGS) -c $<
#.f90.o:
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
#	$(FC) $(COMPFLAGS) -c $<


clean:
	rm -f *.o *.a
