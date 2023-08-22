TARGET = PostProcessCQ
default: $(TARGET)
.SUFFIXES: .f90
include system.make
COMMENT = verstr.f90
ECHOSTR = @echo

DEP_OBJECTS = PostProcessCQ.o read_module.o generic_comms.o input_module.o datatypes_module.o global_module.o numbers_module.o timer_module.o local_module.o dimens_module.o io_module.o timer_stdclocks_module.o species_module.o memory_module.o units.module.o process_module.o output_module.o matrix_data_module.o matrix_module.o basic_types.o block_module.o group_module.o grid_index.module.o primary_module.o maxima_module.o construct_module.o cover_module.o pseudopotential_common.o pseudo_tm_info.o pao_format.o splines_module.o functions_module.o stm_module.o fft_fftw3.o ol_ang_coeff_subs.o
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
