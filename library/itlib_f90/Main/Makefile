#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../itdvrlib_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
MODDIR = ../../Modules
.SUFFIXES: .f90 .mod $(SUFFIXES)

FC = fort
FFLAGS = -c -O5
REL = 

SRCS = \
	dvr_dvd.f90 \
	dvd_err.f90 \
	dvd_dat.f90 \
	trilin.f90 \
	init.f90 \
	invec.f90 \
	hvdvr.f90 \
	hsmall.f90 \
	h1v.f90 \
	h2v.f90 \
	h3v.f90 \
	rdiag.f90 \
	frmres.f90 \
	tstovl.f90 \
	vpert.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f90=.o) 

sources: $(SRCS)
	co -l RCS/*.f90

