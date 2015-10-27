#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(MESA_HOME)/Test_Oleg
include $(TOP_DIR)/Makefile.inc
#

# Typing make make_all will build the full suite of codes
# Typing make make_x with x=bspline or farm will do just that directory
# Same thing for the clean command

make_all: FORCE
	cd $(BSPLINE_MODULES) ; $(MAKE) $(MFLAGS)
	cd $(BSPLINE_SOURCES) ; $(MAKE) $(MFLAGS)
	cd $(FARM) ; $(MAKE) $(MFLAGS)

clean_all: FORCE
	cd $(BSPLINE_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(BSPLINE_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(FARM) ; $(MAKE) $(MFLAGS) clean

make_bspline: FORCE
	cd $(BSPLINE_MODULES) ; $(MAKE) $(MFLAGS)
	cd $(BSPLINE_SOURCES) ; $(MAKE) $(MFLAGS)
clean_bspline: FORCE
	cd $(BSPLINE_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(BSPLINE_SOURCES) ; $(MAKE) $(MFLAGS) clean
make_farm: FORCE
	cd $(FARM) ; $(MAKE) $(MFLAGS)
clean_farm: FORCE
	cd $(FARM) ; $(MAKE) $(MFLAGS) clean

FORCE:

