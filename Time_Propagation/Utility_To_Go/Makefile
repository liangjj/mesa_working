#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the varios flags and libraries needed.
#
.SUFFIXES: .f90 .mod $(SUFFIXES)

MAKE = make
MFLAGS = 
RCS = RCS
#
#
Utility_To_Go_DIRECTORY = \
	Character_Manipulation_Subroutines \
	Common_Los_Alamos_Mathematical_Subroutines \
	General_Utility_Subroutines \
	Mathematical_Subroutines 

all: $(Utility_To_Go_DIRECTORY)

.RECURSIVE: $(Utility_To_Go_DIRECTORY) 
$(Utility_To_Go_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd Character_Manipulation_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd Common_Los_Alamos_Mathematical_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd General_Utility_Subroutines ; $(MAKE) $(MFLAGS) directory
	cd Mathematical_Subroutines ; $(MAKE) $(MFLAGS) directory
clean: FORCE
	cd Character_Manipulation_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd Common_Los_Alamos_Mathematical_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd General_Utility_Subroutines ; $(MAKE) $(MFLAGS) clean
	cd Mathematical_Subroutines ; $(MAKE) $(MFLAGS) clean

FORCE:

