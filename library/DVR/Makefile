# Makefile for making the DVR Library
#
# You should set your environment variables F90, F90FLAGS and MD_LIB to the 
# fortran compiler, the flags you set and the location of the machine dependent 
# fortran.  You can define them here as well but I find that less convenient.
#
#
TOP_DIR = $(MESA_LIB)/DVR
include $(TOP_DIR)/Makefile.inc

all: $(DVR) 

.RECURSIVE: $(DVR)
$(DVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(DVR_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_SOURCES) ; $(MAKE) $(MFLAGS) 
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) 

clean: FORCE
	cd $(DVR_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(DVR_MAIN) ; $(MAKE) $(MFLAGS) clean
FORCE:
