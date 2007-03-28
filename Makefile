# -*- mode: makefile; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
#
# $Id: Makefile,v 1.9.2.1 2006/03/07 07:36:40 drb Exp $
#

#General variables
TARGET = Conquest
default: $(TARGET)
.SUFFIXES: .f .f90
COMMENT = verstr.f90

#Useful variables
DATE = `date +"%b%d%y"`
TARNAME = "CQ_VarNSF_"$(DATE)"_"`date +"%H%M"`".tar"
ECHOSTR = @echo -e
SHELL = /bin/sh

#Include system-dependent variables
include system.make

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

#List of all object files
NODE_OBJECTS = main.o $(MATRIX_OBJS) $(COMM_OBJS) $(ODD_OBJS) $(OVERLAP_OBJS) $(ENERGY_OBJS) $(IONICS_OBJS) $(FORCES_OBJS) $(SETGRID_NEW_OBJS) $(PAO2BLIP_OBJS) $(PS_TM_OBJS)
SRCS = $(NODE_OBJECTS:.o=.f90) basic_types.f90 datatypes.module.f90 matrix_data_module.f90 \
  numbers.module.f90


#Dependency rule
deps.obj.inc: $(SRCS)
	./makemake
	sed /"^mpi.o"/D makemake_deps > deps.obj.inc
	touch $(COMMENT)

#Target
$(TARGET) : $(NODE_OBJECTS) 
	$(FC) $(LINKFLAGS) -o $(TARGET) $(NODE_OBJECTS) $(FDF) $(LIBS) 

#.f90.o:
%.o: %.f90
	$(FC) $(COMPFLAGS) -c $<

.f.o:
	$(FC) $(COMPFLAGS) -c $<

#Dependencies
include deps.obj

$(NODE_OBJECTS): 

initial_read.module.o:initial_read.module.f90
	$(FC) $(COMPFLAGS) -c $<

datestamp.f90: $(COMMENT)
	$(ECHOSTR) "module datestamp\n" > datestamp.f90
	$(ECHOSTR) "  implicit none\n" >> datestamp.f90
	$(ECHOSTR) '  character(len=*), parameter :: datestr="'`date`'"' >> datestamp.f90
	cat $(COMMENT) >> datestamp.f90
	$(ECHOSTR) "\nend module datestamp" >> datestamp.f90
	$(FC) $(COMPFLAGS) -c datestamp.f90

tar:
	tar cvf ../$(TARNAME) *.f *.f90 *.obj Makefile* makemake system.make system/*.make fdf/*.f fdf/Makefile fdf/*.h fdf/README FFT/*.f FFT/Makefile template* utilities/*f90 utilities/Makefile
	gzip ../$(TARNAME)

clean:
	rm -f *.o *.mod *~ libfdf.a libgpfa.a *.d work.pc*
	(cd fdf; make -k clean)
	(cd FFT; make -k clean)
	(cd utilities; make -k clean)

very_clean: 
	make clean; make doc_clean; make pure_clean
doc:
	make xhtml; make html
help:
	$(ECHOSTR) "\nConquest Makefile\n-----------------\n"
	$(ECHOSTR) "\t default (make or make Conquest) is to make Conquest\n"
	$(ECHOSTR) "\t make preConquest makes preConquest\n"
	$(ECHOSTR) "\t make tar makes a .tar.gz file in the directory below\n"
	$(ECHOSTR) "\t make xhtml makes xrefs for ROBODoc html documentation"
	$(ECHOSTR) "\t make html makes ROBODoc html documentation"
	$(ECHOSTR) "\t make doc; netscape Conquest_mi.html gives documentation\n"
	$(ECHOSTR) "\t make <file>_pure removes \"!!\" lines from source file <file>"
	$(ECHOSTR) "\t make pure does for all files (puts into subdir pure)\n"
	$(ECHOSTR) "\t make clean removes all object files and emacs ~ files"
	$(ECHOSTR) "\t make doc_clean cleans up ROBODoc html documentation"
	$(ECHOSTR) "\t make pure_clean removes all _pure files"
	$(ECHOSTR) "\t make very_clean does all above"
list:
	$(ECHOSTR) $(NODE_OBJECTS)
	$(ECHOSTR) $(SOURCES)
	$(ECHOSTR) $(HTMLXREFS)
	$(ECHOSTR) $(HTMLDOCS)

SRC2 = basic_types.f90 datatypes.module.f90 matrix_data_module.f90 \
  numbers.module.f90
include Makefile.Doc
