#   @(#)Makefile	1.4   8/2/91
#
# Makefile for orthlib
#
POLYLIB = ../polylib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)

CPPSRCS = \
	cnvtab.F \
	cnvtpt.F \
	cpoly.F \
	drvply.F \
	genwtf.F \
	genq.F \
	genrwt.F \
	kemel.F \
	lancz.F \
	lgngr.F \
	modab.F \
	orpoly.F \
	plywpt.F \
	recur.F \
	scaprd.F \
	todelta.F \
	toquad.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(POLYLIB)

$(POLYLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(POLYLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F