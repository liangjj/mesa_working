#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(MODULES)
include $(TOP_DIR)/Makefile.inc
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAKE = make
MFLAGS = 
RCS = RCS
LIB_DIRECTORY = \
	General_Modules \
	Utility_Modules \
	Mathematics_Modules 

all: $(LIB_DIRECTORY)

.RECURSIVE: $(LIB_DIRECTORY)

$(LIB_DIRECTORY): FORCE 
	cd $@ ; $(MAKE) $(MFLAGS)



directory: FORCE
	cd General_Modules ; $(MAKE) $(MFLAGS) directory
	cd Utility_Modules ; $(MAKE) $(MFLAGS) directory
	cd Mathematics_Modules ; $(MAKE) $(MFLAGS) directory


clean: FORCE
	cd General_Modules ; $(MAKE) $(MFLAGS) clean
	cd Utility_Modules ; $(MAKE) $(MFLAGS) clean
	cd Mathematics_Modules ; $(MAKE) $(MFLAGS) clean

FORCE:

