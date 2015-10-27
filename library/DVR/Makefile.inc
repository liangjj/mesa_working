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
TOP_DIR = $(MESA_LIB)/DVR
DVR_SOURCES = $(TOP_DIR)/Sources
DVR_MODULES = $(TOP_DIR)/Modules
DVR_MAIN = (TOP_DIR)/Main
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) $(MODULE_LIST) -I$(FEDVR_MODULES) -O4
LDFLAGS = $(F90FLAGS)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<
#
# Basic and Sub Directories
#
DVR = \
	$(DVR_MODULES) \
	$(DVR_SOURCES) \
	$(DVR_MAIN)

DVR_LIB = $(TOP_DIR)/dvr_lib_f90.a
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
