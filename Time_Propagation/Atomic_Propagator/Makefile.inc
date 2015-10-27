#
# Included Makefile for DMRG.
#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_Propagator
#include $(TOP_DIR)/Makefile.sys
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) -I$(MODULES) -I$(ATOMIC_PROP_MODULES) -O4
LDFLAGS = $(F90FLAGS)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<

ATOMIC_PROP_MODULES = $(TOP_DIR)/Modules
ATOMIC_PROP_SOURCES = $(TOP_DIR)/Sources
ATOMIC_PROP_MAIN = $(TOP_DIR)/Main
#
PROPAGATION = \
	    $(ATOMIC_PROP_MODULES) \
	    $(ATOMIC_PROP_SOURCES) \
	    $(ATOMIC_PROP_MAIN)
#
# Libraries
#
ATOMIC_PROP_SOURCE_LIB = $(ATOMIC_PROP_SOURCES)/prop_source_lib_f90.a
ATOMIC_PROP_MOD_LIB = $(ATOMIC_PROP_MODULES)/prop_mod_lib_f90.a
LIB = $(ATOMIC_PROP_MOD_LIB) $(ATOMIC_PROP_SOURCE_LIB) \
	$(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)
#
# Driver
#
DRIVER = $(ATOMIC_PROP_MAIN)/atom_prop
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
