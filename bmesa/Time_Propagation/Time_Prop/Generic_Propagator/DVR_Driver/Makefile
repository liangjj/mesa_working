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
# The Main Directory for the propagator code
#
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
#
# The dvr driver directories consisting of the main, sources and modules.
#
DVR_DRIVER_DIR = $(GENERIC_PROPAGATOR_DIR)/DVR_Driver
DVR_DRIVER_MAIN_DIR = $(DVR_DRIVER_DIR)/DVR_Main
DVR_DRIVER_SOURCES_DIR = $(DVR_DRIVER_DIR)/DVR_Sources
DVR_DRIVER_MOD_DIR = $(DVR_DRIVER_DIR)/DVR_Modules
#
DVR_DRIVER_DIRECTORY = \
	$(DVR_DRIVER_MOD_DIR) \
	$(DVR_DRIVER_SOURCES_DIR) \
	$(DVR_DRIVER_MAIN_DIR) 


all: $(DVR_DRIVER_DIRECTORY) 

dvr_driver: $(DVR_DRIVER_DIRECTORY) 

.RECURSIVE: $(DVR_DRIVER_DIRECTORY)
$(DVR_DRIVER_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(DVR_DRIVER_MOD_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_DRIVER_SOURCES_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(DVR_DRIVER_MAIN_DIR) ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd $(DVR_DRIVER_MOD_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_DRIVER_SOURCES_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_DRIVER_MAIN_DIR) ; $(MAKE) $(MFLAGS) clean

FORCE:

