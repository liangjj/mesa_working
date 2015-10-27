#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Test_Array
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
NOLINK = -c
INCLUDE_LIST = $(MODULE_LIST) -I$(ARRAY_MODULE)
#
.f90.o:
	$(F90) $(NOLINK) $(F90FLAGS) $(INCLUDE_LIST) $<
ARRAY_MAIN = $(TOP_DIR)/Main
ARRAY_MODULE = $(TOP_DIR)/Module
ARRAY_MODULE_LIB = $(TOP_DIR)/Module/array_module_lib.a
#
ARRAY = \
	$(ARRAY_MODULE) \
	$(ARRAY_MAIN)
#
# Libraries
#
LOAD_LIB = $(ARRAY_MODULE_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)
#
# Driver
#
DRIVER = $(ARRAY_MAIN)/driver_x 
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
