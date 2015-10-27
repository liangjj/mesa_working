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
TOP_DIR = $(MESA_HOME)/Time_Propagation/Exponential_Propagator
#include $(TOP_DIR)/Makefile.sys
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) $(MODULE_LIST) -I$(EXPONENTIAL_MODULES) -O4
LDFLAGS = $(F90FLAGS)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<

EXPONENTIAL_MODULES = $(TOP_DIR)/Modules
EXPONENTIAL_SOURCES = $(TOP_DIR)/Sources
EXPONENTIAL_MAIN = $(TOP_DIR)/Main
#
PROPAGATION = \
	    $(EXPONENTIAL_MODULES) \
	    $(EXPONENTIAL_SOURCES) \
	    $(EXPONENTIAL_MAIN)
#
# Libraries
#
EXPONENTIAL_SOURCE_LIB = $(EXPONENTIAL_SOURCES)/exponential_source_lib_f90.a
EXPONENTIAL_MOD_LIB = $(EXPONENTIAL_MODULES)/exponential_mod_lib_f90.a
LIB = $(EXPONENTIAL_MOD_LIB) $(EXPONENTIAL_SOURCE_LIB) \
	$(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)
#
# Driver
#
DRIVER = $(EXPONENTIAL_MAIN)/exp_x
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
