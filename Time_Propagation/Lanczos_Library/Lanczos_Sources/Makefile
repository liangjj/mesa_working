#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Lanczos_Library
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Lanczos_Input.f90

RCS = RCS
CO = co -l

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

%.o : %.f90
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(LANCZOS_MODULES) $<

all: $(LANCZOS_SOURCE_LIB)
$(LANCZOS_SOURCE_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(LANCZOS_SOURCE_LIB) $(?:.f90=.o)
	$(RANLIB) $(LANCZOS_SOURCE_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *~
