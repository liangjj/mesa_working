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
include $(TOP_DIR)/Makefile.sys

#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(DMRG_MODULE_DIR) $<
.f90.o:
	$(F90) $(F90FLAGS) $(MOD_INCLUDE) $<

#
# Directories
#
DMRG_MODULE_DIR  = $(TOP_DIR)/Modules
DMRG_DRIVER_DIR  = $(TOP_DIR)/Main

#
# Driver
#
DMRG_DRIVER      = $(DMRG_DRIVER_DIR)/dmrg

#
# Libraries
#
DMRG_MOD_LIB     = $(DMRG_MODULE_DIR)/dmrg_mod_lib_f90.a

#
# Miscellanea
#
RCS = RCS
CO = co -u
RANLIB = ranlib
