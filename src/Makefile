# -*- mode: makefile; mode: font-lock -*-
#
# $Id: Makefile,v 1.9.2.1 2006/03/07 07:36:40 drb Exp $
#

#General variables
TARGET = Conquest
UTILITIES = tools
default: $(TARGET)
.SUFFIXES: .f .f90
COMMENT = verstr.f90

#Useful variables
DATE = `date +"%b%d%y"`
TARNAME = "CQ_VarNSF_"$(DATE)"_"`date +"%H%M"`".tar"
ECHOSTR = @echo 
SHELL = /bin/sh
NOTIMERS_DIR = altver_notimers
NOSTDTIMERS_DIR = altver_nostdtimers
NOLOCTIMERS_DIR = altver_noloctimers

#Default settings (L.Tong 2012/08/29)
MULT_KERN = default
DIAG_DUMMY = 

#Include system-dependent variables
#Get the system name automatically.
SYSTEM = $(shell hostname -d | cut -d. -f1)
ifneq ($(SYSTEM),)
SYSTEM := .$(SYSTEM)
endif
SYSTEM_PATH = system/system$(SYSTEM).make
include $(SYSTEM_PATH)

#Include lists of object files
include matrix.obj
include comm.obj
include odd.obj
include ionics.obj
include overlap.obj
include energy.obj
include forces.obj
include setgrid_new.obj
include pao2blip.obj
include pseudo_tm.obj
include exx.obj 

#List of all object files
NODE_OBJECTS = main.o datatypes_module.o numbers_module.o datestamp.o                   \
               $(MATRIX_OBJS) $(COMM_OBJS) $(ODD_OBJS) $(OVERLAP_OBJS) $(ENERGY_OBJS)   \
               $(IONICS_OBJS) $(FORCES_OBJS) $(SETGRID_NEW_OBJS) $(PAO2BLIP_OBJS)       \
               $(PS_TM_OBJS)  $(EXX_OBJS)
SRCS = $(NODE_OBJECTS:.o=.f90) basic_types.f90 datatypes_module.f90 matrix_data_module.f90 numbers_module.f90 
# List without datestamp (for dependencies)
NODE_OBJECTS_NODS = main.o datatypes_module.o numbers_module.o                               \
               $(MATRIX_OBJS) $(COMM_OBJS) $(ODD_OBJS) $(OVERLAP_OBJS) $(ENERGY_OBJS)   \
               $(IONICS_OBJS) $(FORCES_OBJS) $(SETGRID_NEW_OBJS) $(PAO2BLIP_OBJS)       \
               $(PS_TM_OBJS)  $(EXX_OBJS)
SRCS_NODS = $(NODE_OBJECTS_NODS:.o=.f90) basic_types.f90 datatypes_module.f90 matrix_data_module.f90 numbers_module.f90 

#Dependency rule
deps.obj.inc: $(SRCS_NODS) $(SYSTEM_PATH)
	touch $(COMMENT)
	$(ECHOSTR) "module datestamp" > datestamp.f90
	$(ECHOSTR) "  implicit none" >> datestamp.f90
	$(ECHOSTR) '  character(len=*), parameter :: datestr="'`date "+%Y/%m/%d at %H:%M %z"`'"' >> datestamp.f90
	sed "s/BBB/"`git symbolic-ref HEAD | cut -b 12-`"/" $(COMMENT) | sed "s/RRR/"`git describe --tags`"/" >> datestamp.f90
	$(ECHOSTR) "end module datestamp" >> datestamp.f90
	./makedeps makedeps.txt $^
	sed /"^mpi.o"/D makedeps.txt > deps.obj.inc

#Target
$(TARGET) : $(NODE_OBJECTS)
	$(FC) $(LINKFLAGS) -o $(TARGET) $(NODE_OBJECTS) $(LIBS)
	cp $(TARGET) ../bin

#.f90.o:
%.o: %.f90
	$(FC) $(COMPFLAGS) -c $<

.f.o:
	$(FC) $(COMPFLAGS) -c $<

#Dependencies
include deps.obj

$(NODE_OBJECTS): 

initial_read_module.o:initial_read_module.f90 datestamp.o
	$(FC) $(COMPFLAGS) -c $<

#datestamp.f90: $(COMMENT)
#	$(ECHOSTR) "module datestamp\n" > datestamp.f90
#	$(ECHOSTR) "  implicit none\n" >> datestamp.f90
#	$(ECHOSTR) '  character(len=*), parameter :: datestr="'`date`'"' >> datestamp.f90
#	cat $(COMMENT) >> datestamp.f90
#	$(ECHOSTR) "\nend module datestamp" >> datestamp.f90
#	$(FC) $(COMPFLAGS) -c datestamp.f90

tar:
	tar cvf ../$(TARNAME) *.f *.f90 *.obj Makefile* makemake $(SYSTEM_PATH) system/*.make FFT/*.f FFT/Makefile template* utilities/*f90 utilities/Makefile
	gzip ../$(TARNAME)

clean:
	rm -f *.o *.mod *~  *.d work.pc* deps.obj.inc
#	(cd FFT; make -k clean)
#	(cd utilities; make -k clean)

#clean_lib:
#	rm -f *.a
#	(cd FFT.FFTE; make -k clean)
#	(cd FFT.GPFA; make -k clean)
#	(cd FFT.ISF;  make -k clean)

very_clean: 
	make clean; make doc_clean; make pure_clean
doc:
	make xhtml; make html
help:
	$(ECHOSTR) "\nConquest Makefile\n-----------------\n"
	$(ECHOSTR) "\t default (make or make Conquest) is to make Conquest\n"
	$(ECHOSTR) "\t make tar makes a .tar.gz file in the directory below\n"
	$(ECHOSTR) "\t make xhtml makes xrefs for ROBODoc html documentation"
	$(ECHOSTR) "\t make html makes ROBODoc html documentation"
	$(ECHOSTR) "\t make clean removes all object files and emacs ~ files"
	$(ECHOSTR) "\t make doc_clean cleans up ROBODoc html documentation"
	$(ECHOSTR) "\t make pure_clean removes all _pure files"
	$(ECHOSTR) "\t make very_clean does all above"
list:
	$(ECHOSTR) $(NODE_OBJECTS)
	$(ECHOSTR) $(SOURCES)
	$(ECHOSTR) $(HTMLXREFS)
	$(ECHOSTR) $(HTMLDOCS)

notimers:
	@if [ ! -d $(NOTIMERS_DIR) ]; \
	then \
	  mkdir $(NOTIMERS_DIR); \
	fi
	@for sourcefile in *.f90; \
	do \
	  grep -vi 'timer_' $${sourcefile}|grep -vi '_timer'|grep -vi 'init_timing_system'|\
	  grep -vi time_report|grep -vi time_threshold|grep -vi tmr_rmv001 > $(NOTIMERS_DIR)/$${sourcefile}; \
	done
	@for objectfile in *.obj; \
	do \
	  grep -vi 'timer_' $${objectfile} > $(NOTIMERS_DIR)/$${objectfile}; \
	done
	@rm -f $(NOTIMERS_DIR)/timer_*.f90
	@cp -pfr FFT Makefile Makefile.Doc utilities *.f makemake $(SYSTEM_PATH) $(NOTIMERS_DIR)/
	@echo "New sources without timers are in $(NOTIMERS_DIR)"

nostdtimers:
	@if [ ! -d $(NOSTDTIMERS_DIR) ]; \
	then \
	  mkdir $(NOSTDTIMERS_DIR); \
	fi
	@for sourcefile in *.f90; \
	do \
	  grep -vi '_stdclocks_' $${sourcefile}|grep -vi 'time_report'|\
	  grep -vi 'tmr_std' > $(NOSTDTIMERS_DIR)/$${sourcefile};\
	done
	@for objectfile in *.obj; \
	do \
	  grep -vi 'timer_std' $${objectfile} > $(NOSTDTIMERS_DIR)/$${objectfile}; \
	done
	@rm -f ${NOSTDTIMERS_DIR}/timer_std*.f90
	@cp -pfr FFT Makefile Makefile.Doc utilities *.f makemake $(SYSTEM_PATH) $(NOSTDTIMERS_DIR)/
	@echo "New sources without standard timers are in $(NOSTDTIMERS_DIR)"

noloctimers:
	@if [ ! -d $(NOLOCTIMERS_DIR) ]; \
	then \
	  mkdir $(NOLOCTIMERS_DIR); \
	fi
	@for sourcefile in *.f90; \
	do \
	  grep -vi 'tmr_l' $${sourcefile} > $(NOLOCTIMERS_DIR)/$${sourcefile}; \
	done
	@cp -pfr FFT Makefile Makefile.Doc utilities *.f *.obj makemake $(SYSTEM_PATH) $(NOLOCTIMERS_DIR)/
	@echo "New sources without local timers are in $(NOLOCTIMERS_DIR)"


SRC2 = basic_types.f90 datatypes.module.f90 matrix_data_module.f90 \
  numbers.module.f90
include Makefile.Doc
