#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#
TOP_DIR = $(LIBRARY)/FEDVR_Lib
#
# System-dependent definitions
#
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)

#FLAGS_AND_MODULES = $(LDFLAGS) $(MODULE_LIST) -I$(FEDVR_MODULES)
NOLINK = -c

.F90.o:
	$(F90) $(NOLINK) $(MODULE_LIST) -I$(FEDVR_MODULES) $(CPPFLAGS) $<
.f90.o:
	$(F90) $(NOLINK) $(MODULE_LIST) -I$(FEDVR_MODULES) $<
.f.o:
	$(F90) $(NOLINK) $(MODULE_LIST) -I$(FEDVR_MODULES) $<
.F.o:
	$(F90) $(NOLINK) $(MODULE_LIST) -I$(FEDVR_MODULES) $(CPPFLAGS) $<

FEDVR_MODULES = $(TOP_DIR)/Modules
#
FEDVR = \
	$(FEDVR_MODULES)
#
# Libraries
#
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
