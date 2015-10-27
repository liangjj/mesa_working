#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(MESA_HOME)/Test_Link
include $(TOP_DIR)/Makefile.inc
#

all: $(TEST)

.RECURSIVE: $(TEST)
$(TEST): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(TEST_MODULES) ; $(MAKE) $(MFLAGS)
	cd $(TEST_MAIN) ; $(MAKE) $(MFLAGS)

clean: FORCE
	cd $(TEST_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(TEST_MAIN) ; $(MAKE) $(MFLAGS) clean

FORCE:

