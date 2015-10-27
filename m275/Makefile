#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
SOURCE = m275
SOURCE_LIB=$(MESA_HOME)/$(SOURCE)/tdse_lib.a
SRCS.f90 = 	\
	$(MESA_HOME)/$(SOURCE)/$(SOURCE).f90
LIB_SRCS.f90 = \
	Data.f90 \
	Derived_Types.f90 \
	Time_Independent_Module.f90 \
	Time_Propagation_Subroutines_Module.f90 \
	Crank_Nicholson_Module.f90 \
	SO_Module.f90 

ALL_INCLUDE = $(INCLUDE_LIST) -I$(MESA_HOME)/$(SOURCE)
include $(MESA_HOME)/Makefile_lib.inc
