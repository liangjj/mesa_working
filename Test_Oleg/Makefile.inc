#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Test_Oleg
#
# Suffix rules
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
#
FARM = $(TOP_DIR)/FARM
BSPLINE = $(TOP_DIR)/BSPLINE
FARM = $(TOP_DIR)/FARM
BSPLINE_MODULES=$(BSPLINE)/Modules
BSPLINE_SOURCES=$(BSPLINE)/Sources
#
#
F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) -I$(BSPLINE_MODULES) -O4
LDFLAGS = $(F90FLAGS)

.f90.o:
	$(F90) $(F90FLAGS) -c $<
.f.o:
	$(F90) $(F90FLAGS) -c $<
#
#
# Driver
#
FARM_X = $(FARM)/farm_x
BSPLINE_MODULE_LIB=$(BSPLINE_MODULES)/bspline_mod_lib.a
BSPLINE_SOURCE_LIB=$(BSPLINE_SOURCES)/bspline_source_lib.a
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
