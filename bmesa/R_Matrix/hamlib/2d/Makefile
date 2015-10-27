#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr hamiltonian routines
#
HAMLIB = ../hamlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)
REL = 

CPPSRCS = \
	drvh2.F \
	h2pac.F \
	ham2d.F \
	setp2.F \
	sym2d.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(HAMLIB)
$(HAMLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(HAMLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F*
