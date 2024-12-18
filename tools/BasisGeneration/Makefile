TARGET = MakeIonFiles
default: $(TARGET)
.SUFFIXES: .f90
#Include system-dependent variables
ifneq ($(SYSTEM),)
$(info System is $(SYSTEM))
SYSTEM_PATH = ../../src/system/system.$(SYSTEM).make
else
SYSTEM_PATH=../../src/system/system.make
endif
ifneq ("$(wildcard $(SYSTEM_PATH))","")
$(info Building using system file $(SYSTEM_PATH))
include $(SYSTEM_PATH)
else
$(info Cannot find system file $(SYSTEM_PATH). Please make one,)
$(info using system/system.example.make as an example, or choose)
$(info an existing file from the system directory using make SYSTEM=label)
$(info to select system/system.label.make)
$(error Compilation aborted.)
endif

COMMENT = verstr.f90
ECHOSTR = @echo 
SHELL = /bin/sh

DEP_OBJECTS = MakeIonFile.o read_module.o generic_comms.o input_module.o datatypes.module.o global_module.o numbers.module.o timer_module.o dimens_local_module.o timer_stdclocks_module.o species_module.o memory_module.o units.module.o pseudo_tm_info.o pseudopotential_common.o pseudo_atom_info_module.o pao_format.o spline_module.o functions_module.o write_module.o schro_module.o mesh_module.o periodic_table_module.o radial_xc_${XC_LIBRARY}_module.o
OBJECTS = $(DEP_OBJECTS) datestamp.o

SRCS = $(DEP_OBJECTS:.o=.f90)
deps.obj: $(SRCS) system.make
	touch $(COMMENT)
	$(ECHOSTR) "module datestamp" > datestamp.f90
	$(ECHOSTR) "  implicit none" >> datestamp.f90
	$(ECHOSTR) '  character(len=*), parameter :: datestr="'`date`'"' >> datestamp.f90
	sed "s/BBB/"`git symbolic-ref HEAD | cut -b 12-`"/" $(COMMENT) | sed "s/RRR/"`git describe --tags`"/" >> datestamp.f90
	$(ECHOSTR) "end module datestamp" >> datestamp.f90
	./makedeps deps.obj datestamp.f90 $^

include deps.obj

$(TARGET): $(OBJECTS)
	$(FC) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)
	cp $(TARGET) ../../bin

%.o: %.f90
	$(FC) $(COMPFLAGS) -c $<

clean:
	rm -f *.o *.mod *~
