#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables F90, F90FLAGS and MD_LIB to the fortran compiler, the flags
# you set and the location of the machine dependentfortran.  You can define them here as well but I
# find that less convenient.
#
#
# This is a comment
TOP_DIR = $(MESA_HOME)/m5000
include $(TOP_DIR)/Makefile.inc

all: $(PROPAGATION) 

.RECURSIVE: $(PROPAGATION)
$(PROPAGATION): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(SIL_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(SIL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(SIL_MAIN) ; $(MAKE) $(MFLAGS) 

clean: FORCE
	cd $(SIL_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(SIL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(SIL_MAIN) ; $(MAKE) $(MFLAGS) clean  
	rm -f *.a *~
FORCE:
