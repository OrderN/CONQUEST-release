#
#  To create ambfor, ambmd: type make
#
#  To create the parallel version of these programs type:
#
#  make clean
#  make mpi
#
#  To run the parallel version of ambfor (fi on 4 processors):
#
# /usr/lib/openmpi/1.3.2-gcc/bin/mpirun -mca btl tcp,self -np 4 ./ambfor.mpi
#
# 

CFLAGS = -c -g
#FC = g77 
FC = gfortran
FFLAGS = -c -g -ffast-math -funroll-loops -O3
LIBS =  -lm
LDR = ${FC} -g

OBJS = angle.o bond.o charge.o conn.o enegrd.o getinp.o param.o \
	tors.o vdwlj.o math.o switch.o wrtout.o pbc.o bsolv.o

all:	ambfor ambmd

# MPI:
#
# If molden does not correctly read the binary files written by ambfor/ambmd:
# This can be caused by use of the g77 compiler for molden and the gfortran
# compiler in ambfor/ambmd (implicitly via mpif77).
# In this case adding the -frecord-marker=4 flag may be used to enforce
# binary compatibility
#

OPENMPI = /usr/lib/openmpi/1.4-gcc
#mpi:    FFLAGS +=  -frecord-marker=4
mpi:	FC = ${OPENMPI}/bin/mpif77
mpi:	CC = ${OPENMPI}/bin/mpicc
mpi:	FFLAGS = -c -g -ffast-math -funroll-loops -O3
mpi:	ambfor.mpi ambmd.mpi

ambfor:	ambfor.o qvdw.o alloc.o optimise.o lbfgs.o $(OBJS)
	$(LDR) -o ambfor ambfor.o qvdw.o alloc.o optimise.o lbfgs.o $(OBJS) $(LIBS)

ambmd:	md.o verlet.o qvdw.o allocmd.o $(OBJS)
	$(LDR) -o ambmd md.o verlet.o qvdw.o allocmd.o $(OBJS) $(LIBS)

ambfor.o:	ambfor.f
	${FC} -c -g -ffast-math -funroll-loops -o ambfor.o ambfor.f

allocmd.o:	alloc.c
	$(CC) $(CFLAGS) -DMD -c alloc.c -o allocmd.o

ambfor.mpi:	ambfor_mpi.f qvdw_mpi.f ambfor_mpi.o qvdw_mpi.o alloc.o optimise.o lbfgs.o $(OBJS)
	${FC} -o ambfor.mpi ambfor_mpi.o qvdw_mpi.o alloc.o optimise.o lbfgs.o $(OBJS) $(LIBS)

ambmd.mpi:	md_mpi.f qvdw_mpi.f md_mpi.o verlet_mpi.o qvdw_mpi.o allocmd.o \
		$(OBJS)
	${FC} -o ambmd.mpi md_mpi.o qvdw_mpi.o allocmd.o verlet_mpi.o $(OBJS) $(LIBS)

ambfor_mpi.f:	ambfor.f
	./tompi ambfor.f

md_mpi.f:	md.f
	./tompi md.f

qvdw_mpi.f:	qvdw.f
	./tompi qvdw.f

verlet_mpi.f:	verlet.f
	./tompi verlet.f

clean:
	\rm *.o
