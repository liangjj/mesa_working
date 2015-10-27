#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the varios flags and libraries needed.
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
#
TOP_DIR=$(DVR)
include $(TOP_DIR)/Makefile.inc
MAKE = make
MFLAGS = 
RCS = RCS
#
# The dvr library subdirectories
#
	DVR_MAIN = $(DVR)/dvrlib_f90/Main
	DVR_SPACE = $(DVR)/dvrlib_f90/Space_f90
	DVR_SHARED = $(DVR)/dvrlib_f90/Shared_f90
	DVR_TIME = $(DVR)/dvrlib_f90/Time_f90
	FD = $(DVR)/fdlib_f90/Main
#
DVR_LIBRARY = \
	$(DVR_TIME) \
	$(DVR_SPACE) \
	$(DVR_SHARED) \
	$(DVR_MAIN) \
	$(FD) 

all: $(DVR_LIBRARY)

.RECURSIVE: $(DVR_LIBRARY) 
$(DVR_LIBRARY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) 
	cd $(FD) ; $(MAKE) $(MFLAGS) 
clean: FORCE
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) clean
	cd $(FD) ; $(MAKE) $(MFLAGS) clean

FORCE:

