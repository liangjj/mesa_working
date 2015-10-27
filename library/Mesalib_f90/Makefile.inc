#
# Including Makefiles must define TOP_DIR appropriately.
#
#
#
# System-dependent definitions
#
TOP_DIR = $(MESA_LIB)/Mesalib_f90
#
# Suffix rules
#
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) -I$(MODULES) -O4
LDFLAGS = $(F90FLAGS)
.f90.o:
	$(F90) $(F90FLAGS) -c $<
#
# Directories
#
MODULES = $(TOP_DIR)/Modules
MATH= $(TOP_DIR)/Mathematical_Subroutines
UTIL = $(TOP_DIR)/General_Utility_Subroutines
#
# Libraries
#
MODULE_LIB = $(MODULES)/mod_lib_f90.a
MATH_LIB = $(MATH)/math_lib_f90.a
UTIL_LIB = $(UTIL)/util_lib_f90.a
#
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
