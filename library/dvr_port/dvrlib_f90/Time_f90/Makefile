#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#

DVR_PORT = $(MESA_LIB)/dvr_port
DVRLIB = $(DVR_PORT)/dvrlib_f90.a
MODDIR = $(DVR_PORT)/Modules

.SUFFIXES: .f90 .mod $(SUFFIXES)

AR = ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O4 -I$(MODDIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	fdf.f90 \
	nrmlze.f90 \
	timpts.f90 \
	tlobato.f90 \
	tphy.f90 \
	tplmat.f90 \
	twpt.f90 

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(DVRLIB)
$(DVRLIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

