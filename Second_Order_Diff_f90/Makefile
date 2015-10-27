#  @(#)Makefile	5.1 11/6/94../Generic_Propagator
#
TOP_DIR = $(MESA_HOME)/Second_Order_Diff_f90
include $(TOP_DIR)/Makefile.inc

all: $(SO)

.RECURSIVE: $(SO)
$(SO): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(SO_MODULES) ; $(MAKE) $(MFLAGS)
	cd $(SO_SOURCES) ; $(MAKE) $(MFLAGS)
	cd $(SO_MAIN) ; $(MAKE) $(MFLAGS)

clean: FORCE
	cd $(SO_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(SO_SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(SO_MAIN) ; $(MAKE) $(MFLAGS) clean

FORCE:
