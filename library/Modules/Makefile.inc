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
#
TOP_DIR = $(MODULES)
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)

NOLINK = -c
F90FLAGS = -L$(MKL_PATH) $(MODULE_LIST) -O3
LDFLAGS = $(F90FLAGS)
.f90.o:
	$(F90) $(NOLINK) $(F90FLAGS) $<

#
# Directories
#
GEN_MOD_LIB = $(TOP_DIR)/gen_mod_lib.a
MATH_MOD_LIB = $(TOP_DIR)/math_lib_mod.a
UTIL_MOD_LIB = $(TOP_DIR)/util_lib_mod.a
UTILITY_LIB = $(TOP_DIR)/mod_lib_f90.a
#
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
