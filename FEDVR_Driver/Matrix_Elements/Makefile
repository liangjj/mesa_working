#@(#)Makefile	1.4   8/2/91
#
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR =$(MESA_HOME)/FEDVR_Driver
include $(TOP_DIR)/Makefile.inc
FEDVR_SRC = FEDVR_Driver.f90

all: $(MATRIX_ELEMENT_DRIVER) 
#MD_LIB=$(MD_LIB_STAT_SEQ)

$(MATRIX_ELEMENT_DRIVER): $(FEDVR_SRC:.f90=.o) $(LIB)
$(MATRIX_ELEMENT_DRIVER): 
#	$(F90) $(LDFLAGS) $(MODULE_LIST) -I$(FEDVR_MODULES) \

#	$(F90) $(LDFLAGS)  \

	$(F90) \
	$(FEDVR_SRC:.f90=.o) \
	$(LIB) $(MD_LIB) -o $(MATRIX_ELEMENT_DRIVER)


files: $(FEDVR_SRC)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
