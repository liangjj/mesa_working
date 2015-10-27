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
TOP_DIR = $(MESA_HOME)/Second_Order_Diff_f90
COULOMB_DIRECTORY = $(MESA_HOME)/Coulomb_Functions_f90
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(SO_MODULES) -I$(COULOMB_MODULES) $<

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
COULOMB_MODULES = $(COULOMB_DIRECTORY)/Coulomb_Modules
SO_MODULES = $(TOP_DIR)/SO_Modules
SO_SOURCES = $(TOP_DIR)/SO_Sources
SO_MAIN = $(TOP_DIR)/SO_Main
#
SO = \
	    $(SO_MODULES) \
	    $(SO_SOURCES) \
	    $(SO_MAIN)
#
# Libraries
MOD_LIB = $(MODULES)/mod_lib_f90.a
POT_LIB = $(POTENTIAL)/pot_lib_f90.a
SO_SOURCE_LIB = $(SO_SOURCES)/so_source_lib_f90.a
SO_MOD_LIB = $(SO_MODULES)/so_mod_lib_f90.a
UTIL_LIB = $(UTILITY)/util_lib_f90.a
DVR_LIB = $(DVR)/dvr_lib_f90.a
FD_LIB = $(DVR)/fd_lib_f90.a
COULOMB_MODULES_LIB = $(COULOMB_MODULES)/coulomb_mod_lib_f90.a
#
# Driver
#
SO_DRIVER = $(SO_MAIN)/so
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
