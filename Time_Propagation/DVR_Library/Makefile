#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the varios flags and libraries needed.
#
.SUFFIXES: .f90 .mod $(SUFFIXES)

MAKE = make
MFLAGS = 
RCS = RCS
#
# The dvr library subdirectories
#
DVR_MAIN = $(DVR_LIB_G_DIR)/Main
DVR_SPACE = $(DVR_LIB_G_DIR)/Space_f90
DVR_SHARED = $(DVR_LIB_G_DIR)/Shared_f90
DVR_TIME = $(DVR_LIB_G_DIR)/Time_f90
FD = $(FD_LIB_G_DIR)/Main
#
DVR_DIRECTORY = \
	$(DVR_TIME) \
	$(DVR_SPACE) \
	$(DVR_SHARED) \
	$(DVR_MAIN) \
	$(FD) 

all: $(DVR_DIRECTORY)

.RECURSIVE: $(DVR_DIRECTORY) 
$(DVR_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) directory
	cd $(FD) ; $(MAKE) $(MFLAGS) directory
clean: FORCE
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) clean
	cd $(FD) ; $(MAKE) $(MFLAGS) clean

FORCE:

