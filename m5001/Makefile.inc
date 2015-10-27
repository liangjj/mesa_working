#
# Included Makefile for DMRG.
#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#

#
# System-dependent definitions
#
include $(MESA_HOME)/Makefile.inc
SOURCE = m5001
TOP_DIR = $(MESA_HOME)/$(SOURCE)
MODULE_DIR = $(TOP_DIR)/.mod_dir
OBJECT_DIR = $(TOP_DIR)/.o_dir
LIB_DIR = $(TOP_DIR)/Library
RUN_DIR = $(TOP_DIR)/Run


SIL_MODULES = $(TOP_DIR)/Modules
SIL_SUBROUTINES = $(TOP_DIR)/Subroutines
SIL_MAIN = $(TOP_DIR)/Main
#
PROPAGATION = \
	    $(SIL_MODULES) \
	    $(SIL_SUBROUTINES) \
	    $(SIL_MAIN)
#
# Libraries
#
SIL_SUBROUTINES_LIB = $(LIB_DIR)/sil_subroutines_lib_f90.a
SIL_MOD_LIB = $(LIB_DIR)/sil_mod_lib_f90.a
LOAD_LIB = $(SIL_MOD_LIB) $(SIL_SUBROUTINES_LIB) \
	  $(FEDVR_MOD_LIB) $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)

