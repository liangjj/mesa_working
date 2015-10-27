#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Legendre
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
MODULE_LIST = -I$(MKL_INCLUDE) -I$(LEGENDRE_MODULES)
#
#     I use the Intel compiler and math libraries and define F90 in 
#     my environment variables. This should be altered for your purposes.
#
F90FLAGS = -L$(MKL_PATH) $(MODULE_LIST) -O4
LDFLAGS = $(F90FLAGS)
#
.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<

LEGENDRE_MAIN = $(TOP_DIR)/Main
LEGENDRE_MODULES = $(TOP_DIR)/Modules
LEGENDRE_MODULES_LIB = $(TOP_DIR)/Modules/legendre_modules_lib.a
#
LEGENDRE = \
	$(LEGENDRE_MODULES) \
	$(LEGENDRE_MAIN)
#
# Libraries
#
LIB = $(LEGENDRE_MODULES_LIB)
#
# Driver
#
DRIVER = $(LEGENDRE_MAIN)/driver_x 
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
