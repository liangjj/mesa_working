#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/CC_Prop
include $(TOP_DIR)/Makefile.inc
SRCS =	CC_Prop_Main.f90 

all: $(PROP_DRIVER)
$(PROP_DRIVER): $(SRCS:.f90=.o) \
                  $(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB) \
                  $(UTIL_LIB) $(MD_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
	$(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB) \
        $(UTIL_LIB) $(MD_LIB) \
	-o $(PROP_DRIVER)     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
