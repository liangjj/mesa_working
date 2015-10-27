#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(LIBRARY)/FEDVR_Lib
include $(TOP_DIR)/Makefile.inc
#

all: $(FEDVR)

.RECURSIVE: $(FEDVR)
$(FEDVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS)

clean: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS) clean

FORCE:

modout: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS) modout
modin: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS) modin
