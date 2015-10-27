#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Coulomb_f90
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Coulomb_Module.f90 \
	Coulfg_Module.f90

RCS = RCS
CO = co -l

.mod.o:

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(COULOMB_MODULES) $<

all: $(COULOMB_MODULES_LIB)
$(COULOMB_MODULES_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(COULOMB_MODULES_LIB) $(?:.f90=.o)
	$(RANLIB) $(COULOMB_MODULES_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *~
