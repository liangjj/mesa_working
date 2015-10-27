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
TOP_DIR = $(MESA_HOME)/Coulomb_f90
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(MODULES) -I$(COULOMB_MODULES) $<
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(COULOMB_MODULES) $<

#
# Basic and Sub Directories
#
MODULES  = $(MESA_LIB)/Modules
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
COULOMB_MODULES = $(TOP_DIR)/Coulomb_Modules
COULOMB_MAIN = $(TOP_DIR)/Coulomb_Main
#
COULOMB = \
	    $(COULOMB_MODULES) \
	    $(COULOMB_MAIN)
#
# Libraries
MOD_LIB = $(MODULES)/mod_lib_f90.a
POT_LIB = $(POTENTIAL)/pot_lib_f90.a
COULOMB_MODULES_LIB = $(COULOMB_MODULES)/coulomb_mod_lib_f90.a
UTIL_LIB = $(UTILITY)/util_lib_f90.a
#
# Driver
#
COULOMB_DRIVER = $(COULOMB_MAIN)/coulomb
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
