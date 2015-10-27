#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#

TOP_DIR = $(MESA_LIB)/Mesalib_f90
include $(TOP_DIR)/Makefile.inc


LIB_DIRECTORY = \
	$(MODULES) \
	$(MATH) \
	$(UTIL)


all: $(LIB_DIRECTORY)

.RECURSIVE: $(LIB_DIRECTORY) 
$(LIB_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) directory
	cd $(MATH) ; $(MAKE) $(MFLAGS) directory
	cd $(UTIL) ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATH) ; $(MAKE) $(MFLAGS) clean
	cd $(UTIL) ; $(MAKE) $(MFLAGS) clean

FORCE:

