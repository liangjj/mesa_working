#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr hamiltonian routines
#
HAMLIB = hamlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

CPPFLAGS = -DECPOINTER
FFLAGS = -c $(CPPFLAGS) 
GET = sccs get
REL = 

CPPSRCS = \
	diagh.F \
	drvh.F \
	rdham.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(HAMLIB)
$(HAMLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(HAMLIB) $(?:.F=.o)

sources: $(CPPSRCS)

	co -l RCS/*.F*
