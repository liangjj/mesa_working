#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
DVR_DRIVER_DIR = $(GENERIC_PROPAGATOR_DIR)/DVR_Driver
DVR_DRIVER_SOURCES_DIR = $(DVR_DRIVER_DIR)/DVR_Sources
DVR_DRIVER_MOD_DIR = $(DVR_DRIVER_DIR)/DVR_Modules
DVR_DRIVER_SOURCES_LIB = $(DVR_DRIVER_SOURCES_DIR)/dvr_driver_sources_lib_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

SRCS = \
	$(DVR_DRIVER_SOURCES_DIR)/pt_reg.f90 \
	$(DVR_DRIVER_SOURCES_DIR)/add_absorb.f90 \
	$(DVR_DRIVER_SOURCES_DIR)/ke_reg_fd.f90 \
	$(DVR_DRIVER_SOURCES_DIR)/ke_reg_dvr.f90 \
	$(DVR_DRIVER_SOURCES_DIR)/modify_diag.f90 \
	$(DVR_DRIVER_SOURCES_DIR)/Input_DVR.f90 \
	$(DVR_DRIVER_SOURCES_DIR)/Space_DVR.f90 

RCS = RCS
CO = co -l

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) -I$(DVR_DRIVER_MOD_DIR) $<

all: $(DVR_DRIVER_SOURCES_LIB)
$(DVR_DRIVER_SOURCES_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(DVR_DRIVER_SOURCES_LIB) $(?:.f90=.o)
	$(RANLIB) $(DVR_DRIVER_SOURCES_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm *.o
