#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Lanczos_Library
include $(TOP_DIR)/Makefile.inc
SRCS =	Lanczos_Main.f90 

all: $(LANCZOS_DRIVER)
$(LANCZOS_DRIVER): $(SRCS:.f90=.o) \
                  $(LANCZOS_MOD_LIB) $(LANCZOS_SOURCE_LIB)   \
                  $(MOD_LIB) $(POT_LIB) \
                  $(UTIL_LIB) $(MD_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
                  $(LANCZOS_MOD_LIB) $(LANCZOS_SOURCE_LIB)   \
        $(MOD_LIB) $(POT_LIB) \
        $(UTIL_LIB) $(MD_LIB) \
	-o $(LANCZOS_DRIVER)     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
