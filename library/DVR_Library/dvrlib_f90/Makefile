# Makefile for making the DVR Library
#
# You should set your environment variables F90, F90FLAGS and MD_LIB to the 
# fortran compiler, the flags you set and the location of the machine dependent 
# fortran.  You can define them here as well but I find that less convenient.
#
#
TOP_DIR = $(DVR)
include $(TOP_DIR)/Makefile.inc

all: $((DVR) 

.RECURSIVE: $(DVR)
$(DVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

DVR_make: FORCE
	cd $(DVR) ; $(MAKE) $(MFLAGS) 

DVR_clean: FORCE
	cd $(DVR) ; $(MAKE) $(MFLAGS) clean

FORCE:
