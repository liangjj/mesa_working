#   @(#)Makefile        1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Coulomb_f90
include $(TOP_DIR)/Makefile.inc
SRCS =  Coulomb_Main.f90

all: $(COULOMB_DRIVER)
$(COULOMB_DRIVER):	$(SRCS:.f90=.o) \
			$(COULOMB_MODULES_LIB) \
			$(MOD_LIB) $(UTIL_LIB) $(MD_LIB)
			$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
			$(COULOMB_MODULES_LIB) $(MOD_LIB) $(UTIL_LIB) $(MD_LIB) \
			-o $(COULOMB_DRIVER)



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~

