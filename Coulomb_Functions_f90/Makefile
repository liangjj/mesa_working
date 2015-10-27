#  @(#)Makefile	5.1 11/6/94../Generic_Propagator
#
#
TOP_DIR = $(MESA_HOME)/Coulomb_Functions_f90
include $(TOP_DIR)/Makefile.inc

all: $(COULOMB) 

.RECURSIVE: $(COULOMB)
$(COULOMB): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

coulomb_make: FORCE
	cd $(COULOMB_MODULES) ; $(MAKE) $(MFLAGS) 
	cd $(COULOMB_MAIN) ; $(MAKE) $(MFLAGS) 

clean: FORCE
	cd $(COULOMB_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(COULOMB_MAIN) ; $(MAKE) $(MFLAGS) clean
FORCE:
