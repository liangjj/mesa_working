#   @(#)Makefile	1.4   8/2/91
#
# Makefile for TIMLIB
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
DVR_LIBRARY =  \
$(GENERIC_PROPAGATOR_DIR)/DVR_Library
DVR_LIB_DIR = $(DVR_LIBRARY)/dvrlib_f90
FD_LIB_DIR = $(DVR_LIBRARY)/fdlib_f90
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
DVR_LIB = $(DVR_LIBRARY)/dvr_lib_f90.a
FD_LIB = $(DVR_LIBRARY)/fd_lib_f90.a

AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FFLAGS = -c -O4

RCS = RCS
CO = co -l

SRCS = \
	band3.f90 \
        band5.f90 \
	band7.f90 \
	cpy_3.f90 \
	fd_input.f90 \
	fd_basis.f90 \
	points.f90 


.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(FD_LIB)
$(FD_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(FD_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o
