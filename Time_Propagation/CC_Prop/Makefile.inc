#
# Included Makefile for DMRG.
#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/04/18 15:23:07 mesa Exp $
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/CC_Prop
#include $(TOP_DIR)/Makefile.sys
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) -I$(MODULES) -I$(PROP_MODULES) -O4
LDFLAGS = $(F90FLAGS)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<

#
# Basic and Sub Directories
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
PROP_SOURCE_LIB = $(PROP_SOURCES)/prop_source_lib_f90.a
PROP_MOD_LIB = $(PROP_MODULES)/prop_mod_lib_f90.a
LIB = $(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)

#
# Driver
#
PROP_DRIVER = $(PROP_MAIN)/ccprop
CLEBSCH_DRIVER = $(PROP_MAIN)/cleb
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
