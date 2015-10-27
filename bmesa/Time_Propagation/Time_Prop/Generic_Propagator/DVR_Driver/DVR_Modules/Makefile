#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
DVR_DRIVER_DIR = $(GENERIC_PROPAGATOR_DIR)/DVR_Driver
DVR_DRIVER_MOD_DIR = $(DVR_DRIVER_DIR)/DVR_Modules
DVR_DRIVER_SOURCES_DIR = $(DVR_DRIVER_DIR)/DVR_Sources
DVR_DRIVER_MOD_LIB = $(DVR_DRIVER_MODULES_DIR)/dvr_driver_modules_lib_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib



LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l
MOD_SRCS = \
	$(DVR_DRIVER_MOD_DIR)/regional_module.f90
.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(DVR_DRIVER_MOD_LIB)
$(DVR_DRIVER_MOD_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_DRIVER_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(DVR_DRIVER_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.mod
