
# NOTE: Dependencies between objects should be explicit


# General variables and targets

.SUFFIXES: .f .f90 .cpp
CQ_UT_PSP2PAO = cq_ut_psp2pao
CQ_UT_HTMTREE = cq_ut_htmtree
CURRENT_WD=$(PWD)
CPP = g++

# Default target

default: all

all: $(CQ_UT_PSP2PAO) $(CQ_UT_HTMTREE)

# Include system dependent variables

include ../system.make


#### LISTS OF OBJECTS

### cq_ut_psp2pao

OBJ_LIST_PSP2PAO =    cq_ut_psp2pao_types.o \
                      cq_ut_psp2pao_input.o \
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
                            cq_ut_psp2pao_input.o \
                            cq_ut_psp2pao_psp.o \
                            cq_ut_psp2pao_numeric.o \
                            cq_ut_psp2pao_hartree.o \
                            cq_ut_psp2pao_xc.o \
                            cq_ut_psp2pao_alloc.o
cq_ut_psp2pao_io.o:         cq_ut_psp2pao_global.o \
                            cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_alloc.o \
                            cq_ut_psp2pao_input.o
cq_ut_psp2pao_main.o:       cq_ut_psp2pao_types.o \
                            cq_ut_psp2pao_input.o \
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
                            cq_ut_psp2pao_alloc.o \
                            cq_ut_psp2pao_input.o
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

### cq_ut_htmtree

OBJ_LIST_HTMTREE = cq_ut_htmtree_main.o \
                   cq_ut_general_cTimerRegister.o \
                   cq_ut_general_cTimerList.o

# Dependencies                   

cq_ut_htmtree_main.o :        cq_ut_general_cTimerRegister.o \
                              cq_ut_general_cTimerList.o
cq_ut_general_cTimerList.o:    cq_ut_general_cTimerRegister.o


# General compilation rules

%.o: %.f90
	$(FC) $(COMPFLAGS) -c $<

.cpp.o:
	$(CPP) $(COMPFLAGS_CPP) -c $<

## TARGETS

# PAO basis generator from pseudopotentials (cq_ut_psp2pao)

$(CQ_UT_PSP2PAO): $(OBJ_LIST_PSP2PAO)
	$(FC) $(LINKFLAGS) -o $(CQ_UT_PSP2PAO) $(OBJ_LIST_PSP2PAO) $(FDF_UT)

# Utility to get timing information in html from a Conquest_out

$(CQ_UT_HTMTREE): $(OBJ_LIST_HTMTREE)
	$(CPP) $(LINKFLAGS) -o $(CQ_UT_HTMTREE) $(OBJ_LIST_HTMTREE)

clean:
#	( cd ../fdf; make clean )
#	rm -f *.o *.mod libfdf.a *~
	rm -f *.o *~ *.mod 
