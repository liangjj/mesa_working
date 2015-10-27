#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(MESA_HOME)/Test_Array
include $(TOP_DIR)/Makefile.inc
#

all: $(ARRAY)

.RECURSIVE: $(ARRAY)
$(ARRAY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(ARRAY_MODULE) ; $(MAKE) $(MFLAGS)
	cd $(ARRAY_MAIN) ; $(MAKE) $(MFLAGS)

clean: FORCE
	cd $(ARRAY_MODULE) ; $(MAKE) $(MFLAGS) clean
	cd $(ARRAY_MAIN) ; $(MAKE) $(MFLAGS) clean

FORCE:

