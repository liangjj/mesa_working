#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Broken_Propagator
include $(TOP_DIR)/Makefile.inc
SRCS =	Atomic_Prop_Main.f90
MD_LIB=$(MD_LIB_DYNAM_SEQ)
all: $(DRIVER)
$(DRIVER): $(SRCS:.f90=.o) \
                  $(ATOMIC_PROP_MOD_LIB) $(ATOMIC_PROP_SOURCE_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) \
                  $(UTILITY_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
	$(ATOMIC_PROP_MOD_LIB) $(ATOMIC_PROP_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) \
        $(UTILITY_LIB) $(MD_LIB) \
	-o $(DRIVER)     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
