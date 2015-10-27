#
# Included Makefile for DMRG.
#
# Including Makefiles must define TOP_DIR appropriately.
#
#
#
# System-dependent definitions
#
#
# Suffix rules
#
TOP_DIR = $(DVR)
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
MODULE_LIST = -I$(MKL_INCLUDE) -I$(GENERAL_MODULES) -I$(MATHEMATICS_MODULES) -I$(UTILITY_MODULES)
F90FLAGS = -L$(MKL_PATH) -I$(MODULE_LIST) -O4
LDFLAGS = $(F90FLAGS)

.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<

#
# Basic and Sub Directories
#
DVR_LIB = $(DVR)/dvr_lib_f90.a
FD_LIB = $(DVR)/fd_lib_f90.a
#
# Miscellaneous
#
MAKE = make
MFLAGS = 
RCS = RCS
CO = co -u
RANLIB = ranlib
