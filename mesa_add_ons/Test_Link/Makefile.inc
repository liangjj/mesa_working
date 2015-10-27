#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#

#
# System-dependent definitions
#
TOP_DIR = $(MESA_HOME)/Test_Link
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
MODULE_LIST = -I$(MKL_INCLUDE) -I$(TEST_MODULES)
F90FLAGS = -L$(MKL_PATH) $(MODULE_LIST) -O4
LDFLAGS = $(F90FLAGS)
#
.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<

TEST_MAIN = $(TOP_DIR)/Main
TEST_MODULES = $(TOP_DIR)/Modules
TEST_MODULES_LIB = $(TOP_DIR)/Modules/test_modules_lib.a
#
TEST = \
	$(TEST_MODULES) \
	$(TEST_MAIN)
#
# Libraries
#
LIB = $(TEST_MODULES_LIB)
#
# Driver
#
DRIVER = $(TEST_MAIN)/driver_x 
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
