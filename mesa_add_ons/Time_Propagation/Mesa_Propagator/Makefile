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

PROPAGATION_DIRECTORY = \
	$(MESA_PROPAGATOR_MODULES) \
	$(MESA_PROPAGATOR_SOURCES) \
	$(MESA_PROPAGATOR_MAIN) 


all: $(PROPAGATION_DIRECTORY) 

prop: $(PROPAGATION_DIRECTORY) 

.RECURSIVE: $(PROPAGATION_DIRECTORY)
$(PROPAGATION_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(MESA_PROPAGATOR_MODULES) ; $(MAKE) $(MFLAGS) directory
	cd $(MESA_PROPAGATOR_SOURCES) ; $(MAKE) $(MFLAGS) directory
	cd $(MESA_PROPAGATOR_MAIN) ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd $(MESA_PROPAGATOR_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(MESA_PROPAGATOR_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(MESA_PROPAGATOR_MAIN) ; $(MAKE) $(MFLAGS) clean

FORCE:

