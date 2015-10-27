#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_FEDVR
include $(TOP_DIR)/Makefile.inc

MOD_SRCS = \
	$(ATOMIC_MODULES)/Atomic_Module.f90

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(ATOMIC_MOD_LIB)
$(ATOMIC_MOD_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(ATOMIC_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(ATOMIC_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
