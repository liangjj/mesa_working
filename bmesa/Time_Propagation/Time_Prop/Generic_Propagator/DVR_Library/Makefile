#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the varios flags and libraries needed.
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
FC = ifort
FFLAGS = -c -O4
MAKE = make
MFLAGS = 
RCS = RCS
#
# The Main Directory
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
DVR_LIBRARY = $(GENERIC_PROPAGATOR_DIR)/DVR_Library
DVR_LIB_DIR = $(DVR_LIBRARY)/dvrlib_f90
FD_LIB_DIR = $(DVR_LIBRARY)/fdlib_f90
#
# The dvr library subdirectories
#
DVR_MAIN_DIR = $(DVR_LIB_DIR)/Main
DVR_SPACE_DIR = $(DVR_LIB_DIR)/Space_f90
DVR_SHARED_DIR = $(DVR_LIB_DIR)/Shared_f90
DVR_TIME_DIR = $(DVR_LIB_DIR)/Time_f90
FD_MAIN_DIR = $(FD_LIB_DIR)/Main
#
DVR_DIRECTORY = \
	$(DVR_TIME_DIR) \
	$(DVR_SPACE_DIR) \
	$(DVR_SHARED_DIR) \
	$(DVR_MAIN_DIR) \
	$(FD_MAIN_DIR) 

all: $(DVR_DIRECTORY)

.RECURSIVE: $(DVR_DIRECTORY) 
$(DVR_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(DVR_SPACE_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_SHARED_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_TIME_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_MAIN_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(FD_MAIN_DIR) ; $(MAKE) $(MFLAGS) directory
clean: FORCE
	cd $(DVR_SPACE_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SHARED_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_TIME_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_MAIN_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(FD_MAIN_DIR) ; $(MAKE) $(MFLAGS) clean

FORCE:

