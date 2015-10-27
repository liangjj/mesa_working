#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
# The Main Directory
MOD_DIR = $(MODULES)
DVR_LIB = $(DVR_LIBRARY_G)/dvr_lib_f90.a
RANLIB = ranlib

RCS = RCS
CO = co -l

SRCS = \
	dvr_input.f90 \
	dvr_basis.f90 \
	read_fourier.f90 \
	read_legendre.f90 \
	read_hermite.f90 \
	read_laguerre.f90 \
	read_grid_parameters.f90 \
	read_theta.f90 \
	ptcal.f90

.f90.o :
	$(F90) $(F90FLAGS) -I$(MOD_DIR) $<

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
