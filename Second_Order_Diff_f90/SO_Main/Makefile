#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Second_Order_Diff_f90
include $(TOP_DIR)/Makefile.inc
SRCS = SO_driver.f90

all: $(SO_DRIVER)
$(SO_DRIVER): $(SRCS:.f90=.o) \
			$(SO_MOD_LIB) \
			$(SO_SOURCE_LIB) \
			$(COULOMB_MODULES_LIB) \
			$(DVR_LIB) \
			$(FD_LIB) \
			$(MOD_LIB) \
			$(POTENTIAL_LIB) \
			$(UTILITY_LIB) \
			$(MD_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
			$(SO_MOD_LIB) \
			$(SO_SOURCE_LIB) \
			$(COULOMB_MODULES_LIB) \
			$(DVR_LIB) \
			$(FD_LIB) \
			$(MOD_LIB) \
			$(POTENTIAL_LIB)  \
			$(UTILITY_LIB) \
			$(MD_LIB) \
	-o $(SO_DRIVER)     



files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
