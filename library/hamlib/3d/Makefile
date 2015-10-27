#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr hamiltonian routines
#
HAMLIB = ../hamlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

CPPFLAGS = -DDECPOINTER
FC = ifort -O4 $(CPPFLAGS)
FFLAGS = -c 
GET = sccs get
REL = 

CPPSRCS = \
	drvh3.F \
	h3pac.F \
	setp3.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(HAMLIB)
$(HAMLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(HAMLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F*


