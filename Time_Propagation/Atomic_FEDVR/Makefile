#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the atomic propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_FEDVR
include $(TOP_DIR)/Makefile.inc
#
all: $(MODULES) $(UTILITY) $(POTENTIAL) $(DVR) $(ATOMIC)
#
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

.RECURSIVE: $(ATOMIC)
$(ATOMIC): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

all_make: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) 
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) 
	cd $(FD) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_MAIN) ; $(MAKE) $(MFLAGS) 

auxiliary_make: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS)
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
	$(MFLAGS) 
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) 
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) 
	cd $(FD) ; $(MAKE) $(MFLAGS) 

main_make: FORCE
	cd $(ATOMIC_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(ATOMIC_MAIN) ; $(MAKE) $(MFLAGS) 

all_clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) clean
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) clean
	cd $(FD) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_MAIN) ; $(MAKE) $(MFLAGS) clean

auxiliary_clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(CHARACTER_MANIPULATION_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(COMMON_LOS_ALAMOS_MATHEMATICAL_SUBROUTINES) ; $(MAKE) \
           $(MFLAGS) clean
	cd $(GENERAL_UTILITY_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATHEMATICAL_SUBROUTINES) ; $(MAKE) $(MFLAGS) clean
	cd $(POTENTIAL) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SHARED) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SPACE) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_TIME) ; $(MAKE) $(MFLAGS) clean
	cd $(FD) ; $(MAKE) $(MFLAGS) clean

main_clean: FORCE
	cd $(ATOMIC_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(ATOMIC_MAIN) ; $(MAKE) $(MFLAGS) clean
FORCE:
