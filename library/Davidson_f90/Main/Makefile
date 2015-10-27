#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_LIB)/Davidson_f90
include $(TOP_DIR)/Makefile.inc
SRCS =	Driver.f90 
#LIB_LIST = $(DAVIDSON_MODULES_LIB) $(DAVIDSON_SHARED_LIB) $(DAVIDSON_SOURCE_LIB) 
LIB_LIST = $(DAVIDSON_MODULES_LIB) $(DAVIDSON_SOURCE_LIB) 
.mod.o:

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(DAVIDSON_MODULES) $<
.f90.mod:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(DAVIDSON_MODULES) $<

all: $(DAVIDSON_DRIVER)
$(DAVIDSON_DRIVER): $(SRCS:.f90=.o) \
			$(LIB_LIST) $(MOD_LIB) $(UTILITY_LIB) $(MD_LIB) 
			$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
			$(LIB_LIST) $(MOD_LIB) $(UTILITY_LIB) $(MD_LIB) \
			-o $(DAVIDSON_DRIVER)     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
