#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Second_Order_Diff_f90
include $(TOP_DIR)/Makefile.inc
MOD_SRCS = \
	regional_module.f90

all: $(SO_MOD_LIB)
$(SO_MOD_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(SO_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(SO_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.mod
