#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(MESA_HOME)/FEDVR_Driver
include $(TOP_DIR)/Makefile.inc
#
all: $(FEDVR)

.RECURSIVE: $(FEDVR)
$(FEDVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS)
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS)

clean: FORCE
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS) clean
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS) clean

make_matrix: FORCE
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS)

clean_matrix: FORCE
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS) clean

make_integrals: FORCE
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS) 

clean_integrals: FORCE
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS) clean

clean_library: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS) clean

FORCE:

