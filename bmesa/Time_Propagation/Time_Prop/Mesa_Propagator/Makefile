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
MESA_PROPAGATOR_DIR = ../Mesa_Propagator
MAIN_DIR = $(MESA_PROPAGATOR_DIR)/Main
PROP_MODULES_DIR = $(MESA_PROPAGATOR_DIR)/Prop_Modules
PROP_SOURCES_DIR = $(MESA_PROPAGATOR_DIR)/Prop_Sources

PROPAGATION_DIRECTORY = \
	$(PROP_MODULES_DIR) \
	$(PROP_SOURCES_DIR) \
	$(MAIN_DIR)


all: $(PROPAGATION_DIRECTORY) 

prop: $(PROPAGATION_DIRECTORY) 

.RECURSIVE: $(PROPAGATION_DIRECTORY)
$(PROPAGATION_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(PROP_MODULES_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(PROP_SOURCES_DIR) ; $(MAKE) $(MFLAGS) directory
	cd $(MAIN_DIR) ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd $(PROP_MODULES_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(PROP_SOURCES_DIR) ; $(MAKE) $(MFLAGS) clean
	cd $(MAIN_DIR) ; $(MAKE) $(MFLAGS) clean

FORCE:

