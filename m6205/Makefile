#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
SOURCE = m6203
SOURCE_LIB=$(MESA_HOME)/$(SOURCE)/grid_lib.a
SRCS.f90 = 	\
	$(MESA_HOME)/$(SOURCE)/$(SOURCE).f90

LIB_SRCS.f90 = \
		Data.f90 \
		Grid_Defined_Types.f90 \
		Gauss_Quadrature.f90 \
		Gauss_Angular_Quadrature.f90 \
		Lebedev_Quadrature.f90 \
		Renormalization.f90 \
		Shell_Info.f90

ALL_INCLUDE = $(INCLUDE_LIST) -I$(MESA_HOME)/$(SOURCE)
include $(MESA_HOME)/Makefile_lib.inc