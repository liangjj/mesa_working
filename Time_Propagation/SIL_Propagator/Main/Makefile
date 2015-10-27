#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/SIL_Propagator
include $(TOP_DIR)/Makefile.inc
SRCS =	SIL_Main.f90
MD_LIB=$(MD_LIB_DYNAM_SEQ)
all: $(DRIVER)
$(DRIVER): $(SRCS:.f90=.o) \
                  $(SIL_MOD_LIB) $(SIL_SOURCE_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) \
                  $(UTILITY_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
	$(SIL_MOD_LIB) $(SIL_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) \
        $(UTILITY_LIB) $(MD_LIB) \
	-o $(DRIVER)     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
