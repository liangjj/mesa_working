#
#
# Including Makefiles must define TOP_DIR appropriately.
#
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_FEDVR
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
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
DVR_LIBRARY = $(MESA_LIB)/DVR_Library
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
DVR = \
	$(DVR_MAIN) \
	$(DVR_SHARED) \
	$(DVR_SPACE) \
	$(DVR_TIME) \
	$(FD) 
#
#
ATOMIC_MODULES = $(TOP_DIR)/Atomic_Modules
ATOMIC_SOURCES = $(TOP_DIR)/Atomic_Sources
ATOMIC_MAIN = $(TOP_DIR)/Atomic_Main
#
#
ATOMIC = \
	$(ATOMIC_MODULES) \
	$(ATOMIC_SOURCES) \
	$(ATOMIC_MAIN)
#
# Libraries
MOD_LIB = $(MODULES)/mod_lib_f90.a
POT_LIB = $(POTENTIAL)/pot_lib_f90.a
UTIL_LIB = $(UTILITY)/util_lib_f90.a
DVR_LIB = $(DVR_LIBRARY)/dvr_lib_f90.a
FD_LIB = $(DVR_LIBRARY)/fd_lib_f90.a
ATOMIC_SOURCE_LIB = $(ATOMIC_SOURCES)/atomic_source_lib_f90.a
ATOMIC_MOD_LIB = $(ATOMIC_MODULES)/atomic_mod_lib_f90.a
#
# Driver
#
ATOMIC_DRIVER = $(ATOMIC_MAIN)/atomic_prop
#
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(ATOMIC_MODULES) $<
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
