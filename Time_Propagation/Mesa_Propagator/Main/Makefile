#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MESA_LIBRARY = $(MESA_LIB)/mesalib.a

RCS = RCS
CO = co -l

SRCS = prop_main.f90

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

%.o : %.f90
	$(F90) $(F90FLAGS) -I$(MOD_DIR) -I$(MESA_PROPAGATOR_MODULES) $<

all: $(MESA_BIN)/m1235
$(MESA_BIN)/m1235: $(SRCS:.f90=.o) \
			$(MESA_PROPAGATOR_LIB) \
			$(DVR_LIB) \
			$(FD_LIB) \
			$(MOD_LIB) \
			$(POTENTIAL_LIB) \
			$(MESA_LIBRARY) \
			$(MD_LIB)
		$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
			$(MESA_PROPAGATOR_LIB) \
			$(DVR_LIB) \
			$(FD_LIB) \
			$(MOD_LIB) \
			$(POTENTIAL_LIB)  \
			$(MESA_LIBRARY) \
			$(MD_LIB) \
		-o $(MESA_BIN)/m1235     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
