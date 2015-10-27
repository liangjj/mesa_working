#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#
#
# System-dependent definitions
#
TOP_DIR = $(UTILITY)
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
NOLINK = -c
#
.F90.o:
	$(F90) $(NOLINK) $(F90FLAGS) $(CPPFLAGS) $<
.f90.o:
	$(F90) $(NOLINK) $(F90FLAGS) $<
.f.o:
	$(F90) $(NOLINK) $(F90FLAGS) $<
.F.o:
	$(F90) $(NOLINK) $(F90FLAGS) $(CPPFLAGS) $<
#
RANLIB = ranlib
FULL_DIR = $(TOP_DIR)/Full_Library
CHRLIB = $(TOP_DIR)/chrlib.a
CLAMLIB = $(TOP_DIR)/clamlib.a
BLASLIB = $(TOP_DIR)/blaslib.a
LAPACKLIB = $(TOP_DIR)/lapacklib.a
DPINTLIB = $(TOP_DIR)/dpintlib.a
INTLIB = $(TOP_DIR)/intlib.a
GENUTILLIB = $(TOP_DIR)/genutil_lib_f90.a
UTILITY_LIB = $(TOP_DIR)/util_lib_f90.a
IOLIB = $(TOP_DIR)/iolib.a
MACHLIB = $(TOP_DIR)/machlib.a
MATHLIB = $(TOP_DIR)/mathlib.a
SYMLIB = $(TOP_DIR)/symlib.a
#
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
