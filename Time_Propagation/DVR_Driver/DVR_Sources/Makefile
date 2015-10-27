#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
DVR_DRIVER_SOURCES = $(DVR_DRIVER)/DVR_Sources
DVR_DRIVER_MOD = $(DVR_DRIVER)/DVR_Modules
DVR_DRIVER_LIB = $(DVR_DRIVER)/dvr_driver_lib_f90.a
RANLIB = ranlib

SRCS = \
	$(DVR_DRIVER_SOURCES)/pt_reg.f90 \
	$(DVR_DRIVER_SOURCES)/add_absorb.f90 \
	$(DVR_DRIVER_SOURCES)/ke_reg_fd.f90 \
	$(DVR_DRIVER_SOURCES)/ke_reg_dvr.f90 \
	$(DVR_DRIVER_SOURCES)/modify_diag.f90 \
	$(DVR_DRIVER_SOURCES)/Input_DVR.f90 \
	$(DVR_DRIVER_SOURCES)/Space_DVR.f90 

RCS = RCS
CO = co -l

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(DVR_DRIVER_MOD) $<

all: $(DVR_DRIVER_LIB)
$(DVR_DRIVER_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(DVR_DRIVER_LIB) $(?:.f90=.o)
	$(RANLIB) $(DVR_DRIVER_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm *.o
