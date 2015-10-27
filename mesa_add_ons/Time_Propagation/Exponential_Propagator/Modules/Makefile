#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Exponential_Propagator
include $(TOP_DIR)/Makefile.inc
MOD_SRCS = \
	Time_Integrator_Derived_Types.f90 \
	Runge_Kutta.f90



all: $(EXPONENTIAL_MOD_LIB)
$(EXPONENTIAL_MOD_LIB): $(MOD_SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(EXPONENTIAL_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(EXPONENTIAL_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
