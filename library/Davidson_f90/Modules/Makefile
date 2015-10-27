#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_LIB)/Davidson_f90
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Davidson_Module.f90



.mod.o:

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) $<
.f90.mod:
	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(DAVIDSON_MODULES_LIB)
$(DAVIDSON_MODULES_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(DAVIDSON_MODULES_LIB) $(?:.f90=.o)
	$(RANLIB) $(DAVIDSON_MODULES_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *~
