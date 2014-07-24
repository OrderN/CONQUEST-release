FC77?=mpif77
FC90?=mpif90
CFLAGS=
FFLAGS?= -fast $(CFLAGS)
LDFLAGS?= $(CFLAGS)
ARFL?=

OBJS =	zfft1d.o zfft3d.o fft235.o kernel.o factor.o

libffte.a: $(OBJS)
	ar $(ARFL) rv libffte.a $(OBJS)

.f.o:	
	$(FC77) $(FFLAGS) -c $<
%.o: %.f90
	$(FC90) $(FFLAGS) -c $<

clean:
	rm -f *.o *.a
