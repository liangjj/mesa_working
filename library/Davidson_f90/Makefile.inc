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
TOP_DIR = $(MESA_LIB)/Davidson_f90
#
# Suffix rules
#
MODULES_LOCATION = -I$(MODULES) -I$(DAVIDSON_MODULES)
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(MODULES) $(MODULES_LOCATION) $<
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(IT_SYS)  $<

DAVIDSON_MODULES = $(TOP_DIR)/Modules
DAVIDSON_SHARED = $(TOP_DIR)/Shared
DAVIDSON_SOURCES = $(TOP_DIR)/Sources
DAVIDSON_MAIN = $(TOP_DIR)/Main
#
DAVIDSON = \
	    $(DAVIDSON_MODULES) \
	    $(DAVIDSON_SHARED) \
	    $(DAVIDSON_SOURCES) \
	    $(DAVAIDSON_MAIN)
#
# Libraries
DAVIDSON_MODULES_LIB = $(DAVIDSON_MODULES)/davidson_modules_lib_f90.a
DAVIDSON_SHARED_LIB = $(DAVIDSON_SHARED)/davidson_shared_lib_f90.a
DAVIDSON_SOURCE_LIB = $(DAVIDSON_SOURCES)/davidson_source_lib_f90.a
#
# Driver
#
DAVIDSON_DRIVER = $(DAVIDSON_MAIN)/davidson
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
