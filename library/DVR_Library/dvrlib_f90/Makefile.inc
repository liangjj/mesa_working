#
# Included Makefile for DMRG.
#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/04/18 15:23:07 mesa Exp $
#

#
# System-dependent definitions
#
TOP_DIR = $(DVR)
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(MODULES) $<
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) $<

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
