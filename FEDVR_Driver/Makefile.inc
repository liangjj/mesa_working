#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#
#
TOP_DIR = $(MESA_HOME)/FEDVR_Driver
# System-dependent definitions
#
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
NOLIST = -c
.F90.o:
	$(F90) $(NOLIST) $(F90FLAGS) $(MODULE_LIST) -I$(FEDVR_MODULES) $(CPPFLAGS) $<
.f90.o:
	$(F90) $(NOLIST) $(F90FLAGS) $(MODULE_LIST) -I$(FEDVR_MODULES) $<

MATRIX_ELEMENTS = $(TOP_DIR)/Matrix_Elements
TWO_ELECTRON_INTEGRALS = $(TOP_DIR)/Two_Electron_Integrals
#
FEDVR = \
	$(MATRIX_ELEMENTS) \
	$(TWO_ELECTRON_INTEGRALS) 
#
# Modules and Libraries
#
FEDVR_MODULES = $(LIBRARY)/FEDVR_Lib/.mod_dir
FEDVR_MODULE_LIB = $(LIBRARY)/FEDVR_Lib/Library/fedvr_mod_lib_f90.a
LIB = $(FEDVR_MODULE_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)
#
# Driver
#
MATRIX_ELEMENT_DRIVER = $(MATRIX_ELEMENTS)/fedvr_me.x 
TWO_ELECTRON_DRIVER = $(TWO_ELECTRON_INTEGRALS)/fedvr_two_el.x 
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
