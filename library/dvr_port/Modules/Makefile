#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#

.SUFFIXES: .f90 .mod $(SUFFIXES)
MODLIB = ../modlib_f90.a
MODDIR = $(MESA_LIB)/dvr_port/Modules

AR = ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O4 -I$(MODDIR)

LD = fort
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	inout.f90 \
	potential.f90 \
	grid_global.f90 \
	prop_global.f90 \
	prop_prnt.f90 \
	dvr_shared.f90 \
	dvr_global.f90 \
	dvr_prnt.f90 \
	dvrprop_global_rt.f90 \
	fd_prnt.f90 \
	fd_global.f90 


.f90.o:
	$(FC) $(FFLAGS) $<

all: $(MODLIB) 
$(MODLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(MODLIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

