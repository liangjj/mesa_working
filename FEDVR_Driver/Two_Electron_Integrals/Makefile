# @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR =$(MESA_HOME)/FEDVR_Driver
include $(TOP_DIR)/Makefile.inc
TWO_ELECTRON_SRC = \
	Two_Electron_Driver.f90 

all: $(TWO_ELECTRON_DRIVER) 
#MD_LIB=$(MD_LIB_STAT_SEQ)

$(TWO_ELECTRON_DRIVER) : $(TWO_ELECTRON_SRC:.f90=.o) $(LIB)
	$(F90) $(LDFLAGS) $(MODULE_LIST) -I$(FEDVR_MODULES) \
	$(TWO_ELECTRON_SRC:.f90=.o) \
	$(LIB) $(MD_LIB) $(GOTOBLASLIB) -o $(TWO_ELECTRON_DRIVER)


files: $(TWO_ELECTRON_SRC)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
