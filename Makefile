TARGET = MakeIonFiles
default: $(TARGET)
.SUFFIXES: .f90
include system.make

OBJECTS = MakeIonFile.o read_module.o generic_comms.o input_module.o datatypes.module.o global_module.o numbers.module.o timer_module.o dimens_local_module.o timer_stdclocks_module.o species_module.o memory_module.o units.module.o pseudo_tm_info.o pseudopotential_common.o pao_info_module.o pao_format.o spline_module.o functions_module.o write_module.o schro_module.o mesh_module.o periodic_table_module.o radial_xc_module.o

SRCS = $(OBJECTS:.o=.f90)
deps.obj: $(SRCS) system.make
	./makedeps deps.obj $^

include deps.obj

$(TARGET): $(OBJECTS)
	$(FC) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

%.o: %.f90
	$(FC) $(COMPFLAGS) -c $<
