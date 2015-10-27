#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
DVR_DRIVER_MOD = $(DVR_DRIVER)/DVR_Modules
DVR_DRIVER_SOURCES = $(DVR_DRIVER)/DVR_Sources
DVR_DRIVER_LIB = $(DVR_DRIVER)/dvr_driver_lib_f90.a
RANLIB = ranlib

RCS = RCS
CO = co -l
MOD_SRCS = \
	$(DVR_DRIVER_MOD)/regional_module.f90
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(DVR_DRIVER_LIB)
$(DVR_DRIVER_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_DRIVER_LIB) $(?:.f90=.o)
	$(RANLIB) $(DVR_DRIVER_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.mod
