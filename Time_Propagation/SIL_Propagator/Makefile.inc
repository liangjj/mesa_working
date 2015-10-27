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
TOP_DIR = $(MESA_HOME)/Time_Propagation/SIL_Propagator
#include $(TOP_DIR)/Makefile.sys
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) $(MODULE_LIST) -I$(SIL_MODULES) -O3
LDFLAGS = $(F90FLAGS)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<

SIL_MODULES = $(TOP_DIR)/Modules
SIL_SOURCES = $(TOP_DIR)/Sources
SIL_MAIN = $(TOP_DIR)/Main
#
PROPAGATION = \
	    $(SIL_MODULES) \
	    $(SIL_SOURCES) \
	    $(SIL_MAIN)
#
# Libraries
#
SIL_SOURCE_LIB = $(SIL_SOURCES)/sil_source_lib_f90.a
SIL_MOD_LIB = $(SIL_MODULES)/sil_mod_lib_f90.a
LIB = $(SIL_MOD_LIB) $(SIL_SOURCE_LIB) \
	$(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)
#
# Driver
#
DRIVER = $(SIL_MAIN)/sil_x
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
