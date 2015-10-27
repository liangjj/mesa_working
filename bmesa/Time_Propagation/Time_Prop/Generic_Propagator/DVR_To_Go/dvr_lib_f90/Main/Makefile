#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
DVR_LIB = ../../dvr_lib_f90.a
PROP_MOD_DIR = ../../../Prop_To_Go/Modules
MOD_DIR = ../../../Modules
AR = ar
ARFLAGS = rvu
RANLIB = ranlib

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

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.a

makercs:
	ci -l *.f90
