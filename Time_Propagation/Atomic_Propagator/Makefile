#
# Makefile for making the propagator codes
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables F90, F90FLAGS and MD_LIB to the fortran compiler, the flags
# you set and the location of the machine dependentfortran.  You can define them here as well but I
# find that less convenient.
#
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_Propagator
include $(TOP_DIR)/Makefile.inc

all: $(MODULES) $(UTILITY) $(POTENTIAL) $(DVR) $(PROPAGATION) 

.RECURSIVE: $(MODULES)
$(MODULES): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(UTILITY)
$(UTILITY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(POTENTIAL)
$(POTENTIAL): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(DVR)
$(DVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

.RECURSIVE: $(PROPAGATION)
$(PROPAGATION): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

all_make: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) 
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) 
	cd $(DVR) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_PROP_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_PROP_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_PROP_MAIN) ; $(MAKE) $(MFLAGS) 

module_make: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 
libraries_make: FORCE
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) 
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) 
DVR_make: FORCE
	cd $(DVR) ; $(MAKE) $(MFLAGS) 

main_make: FORCE
	cd $(ATOMIC_PROP_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_PROP_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_PROP_MAIN) ; $(MAKE) $(MFLAGS) 

all_clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) clean
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_PROP_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_PROP_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_PROP_MAIN) ; $(MAKE) $(MFLAGS) clean

module_clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
libraries_clean: FORCE
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) clean
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) clean
DVR_clean: FORCE
	cd $(DVR) ; $(MAKE) $(MFLAGS) clean

main_clean: FORCE
	cd $(ATOMIC_PROP_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_PROP_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_PROP_MAIN) ; $(MAKE) $(MFLAGS) clean  

FORCE:
