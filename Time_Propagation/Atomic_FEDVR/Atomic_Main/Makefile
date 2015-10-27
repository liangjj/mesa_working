#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the atomic propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_FEDVR
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Atomic_Main.f90 \

.f90.0:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(ATOMIC_MODULES) $<

all: $(BIN)/Atomic
$(BIN)/Atomic: $(SRCS:.f90=.o) \
                  $(ATOMIC_MOD_LIB) (ATOMIC_SOURCE_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB) \
                  $(UTIL_LIB) $(MD_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
        $(ATOMIC_MOD_LIB) $(ATOMIC_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB)  \
	$(UTIL_LIB) $(MD_LIB) \
        -o $(BIN)/Atomic     
