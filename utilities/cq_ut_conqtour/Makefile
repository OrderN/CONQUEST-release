.PHONY: manual clean deepclean deepdeepclean

OBJECTS=cq_ut_conqtour_main.o\
        cq_ut_conqtour_init.o\
        cq_ut_conqtour_control.o\
        cq_ut_conqtour_io.o\
        cq_ut_conqtour_keyboard.o\
        cq_ut_conqtour_mouse.o\
        cq_ut_conqtour_legend.o\
        cq_ut_conqtour_glview.o\
        cq_ut_conqtour_img.o\
        cq_ut_conqtour_pymol.o\
        cq_ut_conqtour_slicevis.o\
        cq_ut_conqtour_slicelog.o\
        cq_ut_conqtour_math.o\
        cq_ut_conqtour_glthread.o\
        cq_ut_conqtour_pythread.o\
        cq_ut_conqtour_pymodule.o
MAINTARGET=cq_ut_conqtour
HEADER=cq_ut_conqtour.h

ifdef CQ_SYSTEM
  include system.make.$(CQ_SYSTEM)
else
  include system.make
endif

%.o: %.c $(HEADER)
	$(CC) $(USRFLAGS) $(CFLAGS) $(INCDIRS) -c $<

default: $(MAINTARGET) manual

$(MAINTARGET): $(OBJECTS)
	$(CC) $(USRFLAGS) $(CFLAGS) $(LIBDIRS) $(LIBS) $(OBJECTS) -o $(MAINTARGET)
ifdef CQ_USRBIN
	cp -p $(MAINTARGET) $(CQ_USRBIN)/$(MAINTARGET)
endif

manual:
	$(MAKE) -C manual
	
clean:
	rm -f $(OBJECTS)
	$(MAKE) -C manual clean
	
deepclean:
	rm -f $(OBJECTS) $(MAINTARGET)
	$(MAKE) -C manual deepclean

deepdeepclean:
	rm -f $(OBJECTS) $(MAINTARGET) $(USRBIN)/$(MAINTARGET)
