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
# The Main Directory for the propagator code
#
# The dvr driver directories consisting of the main, sources and modules.
#
DVR_DRIVER = $(B_MESA)/DVR_Driver
DVR_DRIVER_MAIN = $(DVR_DRIVER)/DVR_Main
DVR_DRIVER_SOURCES = $(DVR_DRIVER)/DVR_Sources
DVR_DRIVER_MOD = $(DVR_DRIVER)/DVR_Modules
#
DVR_DRIVER_DIRECTORY = \
	$(DVR_DRIVER_MOD) \
	$(DVR_DRIVER_SOURCES) \
	$(DVR_DRIVER_MAIN) 


all: $(DVR_DRIVER_DIRECTORY) 

dvr_driver: $(DVR_DRIVER_DIRECTORY) 

.RECURSIVE: $(DVR_DRIVER_DIRECTORY)
$(DVR_DRIVER_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(DVR_DRIVER_MOD) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_DRIVER_SOURCES) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_DRIVER_MAIN) ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd $(DVR_DRIVER_MOD) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_DRIVER_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_DRIVER_MAIN) ; $(MAKE) $(MFLAGS) clean

FORCE:

