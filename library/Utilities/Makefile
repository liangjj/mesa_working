#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#

#
LIB_DIRECTORY = \
	Character_Manipulation_Subroutines \
	IOsys_Subroutines \
	Mathematical_Subroutines \
	Common_Los_Alamos_Mathematical_Subroutines \
	General_Utility_Subroutines \
	Machine_Dependent_Subroutines/opteron \
	Integral_Library \
	Double_Precision_Integral_Library \
	Symmetry_Library

all: $(LIB_DIRECTORY)

.RECURSIVE: $(LIB_DIRECTORY) 
$(LIB_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd Character_Manipulation_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd IOsys_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd Mathematical_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd Common_Los_Alamos_Mathematical_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd Machine_Dependent_Subroutines/opteron ; $(MAKE) $(MFLAGS) directory
	cd General_Utility_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd Integral_Library ; $(MAKE) $(MFLAGS) directory
	cd Double_Precision_Integral_Library ; $(MAKE) $(MFLAGS) directory
	cd Symmetry_Library ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd Character_Manipulation_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd IOsys_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd Mathematical_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd Common_Los_Alamos_Mathematical_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd Machine_Dependent_Subroutines/opteron ; $(MAKE) $(MFLAGS) clean
	cd General_Utility_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd Integral_Library ; $(MAKE) $(MFLAGS) clean
	cd Double_Precision_Integral_Library ; $(MAKE) $(MFLAGS) clean
	cd Symmetry_Library ; $(MAKE) $(MFLAGS) clean

FORCE:

