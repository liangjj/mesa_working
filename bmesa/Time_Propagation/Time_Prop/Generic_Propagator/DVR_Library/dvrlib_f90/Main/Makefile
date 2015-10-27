#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
# The Main Directory
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
DVR_LIBRARY =  \
$(GENERIC_PROPAGATOR_DIR)/DVR_Library
DVR_LIB_DIR = $(DVR_LIBRARY)/dvrlib_f90
FD_LIB_DIR = $(DVR_LIBRARY)/fdlib_f90
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
DVR_LIB = $(DVR_LIBRARY)/dvr_lib_f90.a
FD_LIB = $(DVR_LIBRARY)/fd_lib_f90.a

AR = ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O4 -I$(MOD_DIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	atomic_input.f90 \
	atomic_basis.f90 \
	dvr_input.f90 \
	dvr_basis.f90 \
	read_fourier.f90 \
	read_legendre.f90 \
	read_grid_parameters.f90 \
	ptcal.f90

.f90.o :
	$(FC) $(FFLAGS) $<

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

