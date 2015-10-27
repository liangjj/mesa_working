#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
SOURCE = m6204
GRID_LIB_DIR=$(LIBRARY)/gridlib
GRID_LIB=$(GRID_LIB_DIR)/grid_lib.a
#
#   Force the library to be made if needed.
#

all: $(GRID_LIB_DIR)
$(GRID_LIB_DIR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)
FORCE:

SOURCE_LIB=$(GRID_LIB)
SRCS.f90 = 	\
	$(MESA_HOME)/$(SOURCE)/$(SOURCE).f90

ALL_INCLUDE = $(INCLUDE_LIST) -I$(MESA_HOME)/$(SOURCE) -I$(GRID_LIB_DIR)
include $(MESA_HOME)/Makefile_lib.inc