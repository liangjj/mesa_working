#
# Included Makefile for DMRG.
#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id$
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Generic_Propagator
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(MODULES) -I$(PROP_MODULES) $<
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(PROP_MODULES) $<

#
# Basic and Sub Directories
#
MODULES  = $(MESA_LIB)/Modules
POTENTIAL = $(MESA_LIB)/Potential
UTILITY = $(MESA_LIB)/Utilities
	CHARACTER_MANIPULATION_SUBROUTINES = \
		$(UTILITY)/Character_Manipulation_Subroutines 
	COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES = \
		$(UTILITY)/Common_Los_Alamos_Mathematical_Subroutines 
	GENERAL_UTILITY_SUBROUTINES = \
		$(UTILITY)/General_Utility_Subroutines 
	MATHEMATICAL_SUBROUTINES = \
		$(UTILITY)/Mathematical_Subroutines
DVR = $(MESA_LIB)/DVR_Library
	DVR_MAIN = \
		$(DVR_LIBRARY)/dvrlib_f90/Main 
	DVR_SHARED = \
		$(DVR_LIBRARY)/dvrlib_f90/Shared_f90
	DVR_SPACE = \
		$(DVR_LIBRARY)/dvrlib_f90/Space_f90 
	DVR_TIME = \
		$(DVR_LIBRARY)/dvrlib_f90/Time_f90
	FD = \
		$(DVR_LIBRARY)/fdlib_f90/Main
#
#
PROP_MODULES = $(TOP_DIR)/Modules
PROP_SOURCES = $(TOP_DIR)/Sources
PROP_MAIN = $(TOP_DIR)/Main
#
PROPAGATION = \
	    $(PROP_MODULES) \
	    $(PROP_SOURCES) \
	    $(PROP_MAIN)
#
# Libraries
MOD_LIB = $(MODULES)/mod_lib_f90.a
POT_LIB = $(POTENTIAL)/pot_lib_f90.a
PROP_SOURCE_LIB = $(PROP_SOURCES)/prop_source_lib_f90.a
PROP_MOD_LIB = $(PROP_MODULES)/prop_mod_lib_f90.a
UTIL_LIB = $(UTILITY)/util_lib_f90.a
DVR_LIB = $(DVR)/dvr_lib_f90.a
FD_LIB = $(DVR)/fd_lib_f90.a
#
# Driver
#
PROP_DRIVER = $(PROP_MAIN)/dvrprop
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
