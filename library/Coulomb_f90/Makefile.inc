#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#
#
# System-dependent definitions
#
TOP_DIR = ./
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) $(MODULE_LIST) -O3
LDFLAGS = $(F90FLAGS)
#
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.f.o:
	$(F90) $(F90FLAGS) -c $<
.F.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
#
RANLIB = ranlib

COUL_LIB = $(TOP_DIR)/coullib.a
#
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
