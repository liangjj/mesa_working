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
TOP_DIR = $(MESA_HOME)/Time_Propagation/Lanczos_Library
#include $(TOP_DIR)/Makefile.sys
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(MODULES) -I$(LANCZOS_MODULES) $<
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(LANCZOS_MODULES) $<

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
#
#
LANCZOS_MODULES = $(TOP_DIR)/Lanczos_Modules
LANCZOS_SOURCES = $(TOP_DIR)/Lanczos_Sources
LANCZOS_MAIN = $(TOP_DIR)/Main
#
LANCZOS = \
	    $(LANCZOS_MODULES) \
	    $(LANCZOS_SOURCES) \
	    $(LANCZOS_MAIN)
#
# Libraries
MOD_LIB = $(MODULES)/mod_lib_f90.a
POT_LIB = $(POTENTIAL)/pot_lib_f90.a
LANCZOS_SOURCE_LIB = $(LANCZOS_SOURCES)/lanczos_source_lib_f90.a
LANCZOS_MOD_LIB = $(LANCZOS_MODULES)/lanczos_mod_lib_f90.a
UTIL_LIB = $(UTILITY)/util_lib_f90.a
#
# Driver
#
LANCZOS_DRIVER = $(LANCZOS_MAIN)/lanczos
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
