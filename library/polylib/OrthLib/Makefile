#   @(#)Makefile	1.4   8/2/91
#
# Makefile for orthlib
#
ORTHLIB = ../orthlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)

CPPSRCS = \
	drvply.F \
	orpoly.F \
	poly.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(ORTHLIB)

$(ORTHLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(ORTHLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F