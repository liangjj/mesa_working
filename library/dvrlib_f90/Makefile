#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the varios flags and libraries needed.
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
#FC = ifort
#FFLAGS = -c -O4
MAKE = make
MFLAGS = 
RCS = RCS
#
# The Main Directories
DVR_LIBRARY = $(MESA_LIB)/DVR_Library
DVR_LIB_DIR = $(DVR_LIBRARY)/dvrlib_f90
FD_LIB_DIR = $(DVR_LIBRARY)/fdlib_f90
#
# The sub directories
#
DVR_MAIN_DIR = $(DVR_LIB_DIR)/Main
DVR_SPACE_DIR = $(DVR_LIB_DIR)/Space_f90
DVR_SHARED_DIR = $(DVR_LIB_DIR)/Shared_f90
DVR_TIME_DIR = $(DVR_LIB_DIR)/Time_f90
FD_MAIN_DIR = $(FD_LIB_DIR)/Main
#
DVR_DIRECTORY = \
	$(DVR_MAIN_DIR) \
	$(DVR_SPACE_DIR) \
	$(DVR_SHARED_DIR) \
	$(DVR_TIME_DIR) 

FD_DIRECTORY = \
	$(FD_MAIN_DIR) 


all: $(DVR_DIRECTORY) $(FD_DIRECTORY)

dvr: $(DVR_DIRECTORY) 

.RECURSIVE: $(DVR_DIRECTORY) $(FD_DIRECTORY)
$(DVR_DIRECTORY) $(FD_DIRECTORY): FORCE
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

