#  @(#)Makefile	5.1 11/6/94../Generic_Propagator
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables F90, F90FLAGS and MD_LIB to the fortran compiler, the flags
# you set and the location of the machine dependentfortran.  You can define them here as well but I
# find that less convenient.
#
# There is a grand makefile in the CC_Prop directory that with the command make will compile and load
# everything.  Similarly a make clean will remove all but the sources.  Ypu can also go into the individual
# subdirectoies and do individual makes.
#
# The CC_Prop directory contains the sources, modules and main code.  The DVR_Library directory builds the
# FEDVR library which computes everthing needed to compute the one-body operators in the FEDVR basis as well
# as some other things like the Poisson solver for the two-electron integrals.  The Utilities directory is a hack
# of the old Mesalib directories.  Its mostly all of it except  IOsys.  I like to use these routines over and over again
# as they do lots of good things and rewriting code is not my favorite pastime.  This directory would basically allow
# a user to function independently of BLAS, LAPACK etc. if they needed to since it has fortran versions of all of
# that and a heck of a lot more such as the LANL CLAMS library and the SLATEC library.
#
# This Makefile defines the major directories and then recursively goes into each of them in the proper order
# and compiles and loads what is needed.  The last directory, Main, is where the executable is created and where all the
# externals get satisfied.
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/CC_Prop
include $(TOP_DIR)/Makefile.inc

all: $(MODULES) $(UTILITY) $(POTENTIAL) $((DVR) $(PROPAGATION) 

.RECURSIVE: $(MODULES)
$(MODULES): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(UTILITY)
$(UTILITY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(POTENTIAL)
$(POTENTIAL): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(DVR)
$(DVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(PROPAGATION)
$(PROPAGATION): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

all_make: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) 
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) 
	cd $(DVR) ; $(MAKE) $(MFLAGS) 
	cd $(PROP_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(PROP_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(PROP_MAIN) ; $(MAKE) $(MFLAGS) 

module_make: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 
libraries_make: FORCE
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) 
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) 
DVR_make: FORCE
	cd $(DVR) ; $(MAKE) $(MFLAGS) 

main_make: FORCE
	cd $(PROP_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(PROP_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(PROP_MAIN) ; $(MAKE) $(MFLAGS) 

all_clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) clean
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR) ; $(MAKE) $(MFLAGS) clean
	cd $(PROP_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(PROP_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(PROP_MAIN) ; $(MAKE) $(MFLAGS) clean

module_clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
libraries_clean: FORCE
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) clean
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) clean
DVR_clean: FORCE
	cd $(DVR) ; $(MAKE) $(MFLAGS) clean

main_clean: FORCE
	cd $(PROP_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(PROP_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(PROP_MAIN) ; $(MAKE) $(MFLAGS) clean  

FORCE:
