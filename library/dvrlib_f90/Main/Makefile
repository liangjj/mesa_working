#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
DVRLIB = ../dvrlib_f90.a
MODDIR = ../../Modules

.SUFFIXES: .f90 .mod $(SUFFIXES)

AR = ar
ARFLAGS = rvu
RANLIB = ranlib

#FFLAGS = -c -O4 -I$(MODDIR)

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
	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(DVRLIB)
$(DVRLIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.f90=.o)
	$(RANLIB) $(DVRLIB)
sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

