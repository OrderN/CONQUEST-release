
# NOTE: Dependencies between objects should be explicit


# General variables and targets

.SUFFIXES: .f .f90
CQ_UT_PSP2PAO = cq_ut_psp2pao
FDF=libfdf.a
CURRENT_WD=$(PWD)


# Default target

default: $(CQ_UT_PSP2PAO)


# Include system dependent variables

include ../system.make


## LISTS OF OBJECTS

# cq_ut_psp2pao

OBJ_LIST_PSP2PAO = $(FDF) \
                      cq_ut_psp2pao_types.o \
                      cq_ut_psp2pao_global.o \
                      cq_ut_psp2pao_alloc.o \
                      cq_ut_psp2pao_numeric.o \
                      cq_ut_psp2pao_psp.o \
                      cq_ut_psp2pao_density.o \
                      cq_ut_psp2pao_hartree.o \
                      cq_ut_psp2pao_xc.o \
                      cq_ut_psp2pao_initialise.o \
                      cq_ut_psp2pao_io.o \
                      cq_ut_psp2pao_ev.o \
                      cq_ut_psp2pao_solve.o \
                      cq_ut_psp2pao_scf.o \
                      cq_ut_psp2pao_basis.o \
                      cq_ut_psp2pao_orbital.o \
                      cq_ut_psp2pao_main.o

# Dependencies

cq_ut_psp2pao_alloc.o:      cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o
cq_ut_psp2pao_basis.o:      cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_psp.o \
                            cq_ut_psp2pao_numeric.o \
                            cq_ut_psp2pao_hartree.o \
                            cq_ut_psp2pao_xc.o \
                            cq_ut_psp2pao_alloc.o
cq_ut_psp2pao_density.o:    cq_ut_psp2pao_global.o
cq_ut_psp2pao_ev.o:         cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_numeric.o \
                            cq_ut_psp2pao_xc.o \
                            cq_ut_psp2pao_hartree.o \
                            cq_ut_psp2pao_density.o
cq_ut_psp2pao_global.o:     cq_ut_psp2pao_types.o
cq_ut_psp2pao_hartree.o:    cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_numeric.o 
cq_ut_psp2pao_initialise.o: cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_psp.o \
                            cq_ut_psp2pao_numeric.o \
                            cq_ut_psp2pao_hartree.o \
                            cq_ut_psp2pao_xc.o \
                            cq_ut_psp2pao_alloc.o
cq_ut_psp2pao_io.o:         cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_alloc.o
cq_ut_psp2pao_main.o:       $(FDF) \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_initialise.o \
                            cq_ut_psp2pao_scf.o \
                            cq_ut_psp2pao_density.o \
                            cq_ut_psp2pao_io.o \
                            cq_ut_psp2pao_basis.o \
                            cq_ut_psp2pao_orbital.o \
                            cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_alloc.o
cq_ut_psp2pao_numeric.o:    cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_global.o
cq_ut_psp2pao_orbital.o:    cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_numeric.o
cq_ut_psp2pao_psp.o:        cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_numeric.o \
                            cq_ut_psp2pao_alloc.o
cq_ut_psp2pao_scf.o:        cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_solve.o \
                            cq_ut_psp2pao_ev.o
cq_ut_psp2pao_solve.o:      cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_numeric.o
cq_ut_psp2pao_xc.o:         cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_numeric.o

# General compilation rules

%.o: %.f90
	$(FC) $(COMPFLAGS) -c $<

## TARGETS

# FDF Module

fdf:
	@echo $(CURRENT_WD)
	(cd ../fdf; make; cp *.mod $(CURRENT_WD); cp libfdf.a $(CURRENT_WD))

#	@echo $(CURRENT_WD)
#	(cd ../fdf; $(MAKE) "FC=$(FC)" "FFLAGS=$(COMPFLAGS)" module; cp *.mod $(CURRENT_WD); cp libfdf.a $(CURRENT_WD))

# PAO basis generator from pseudopotentials (cq_ut_psp2pao)

$(CQ_UT_PSP2PAO): $(OBJ_LIST_PSP2PAO)
	$(FC) $(LINKFLAGS) -o $(CQ_UT_PSP2PAO) $(OBJ_LIST_PSP2PAO) $(FDF)

clean:
#	( cd ../fdf; make clean )
#	rm -f *.o *.mod libfdf.a *~
	rm -f *.o *~
