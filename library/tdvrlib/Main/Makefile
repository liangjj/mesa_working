#   @(#)Makefile	1.4   8/2/91
#
# Makefile for tdvrlib
#
TDVRLIB = ../tdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)
GET = sccs get

CPPSRCS = \
	addd.F \
	anlyze.F \
	cmpslv.F \
	eigslv.F \
	fil1.F \
	fil2.F \
	fil3.F \
	fildr.F \
	filvt.F \
	gpaket.F \
	hamt.F \
	hamtot.F \
	initpsi.F \
	lindat.F \
	mkphi.F \
	mkpsi0.F \
	mkhpsi.F \
	nospac.F \
	nrpaket.F \
	plotit.F \
	plot1.F \
	pltwfn.F \
	radpakt.F \
	set.F \
	sppose.F \
	statvc.F \
	vcoef1.F \
	vcoef2.F \
	vcoef3.F \
	vect0.F \
	zproj.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(TDVRLIB)
$(TDVRLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(TDVRLIB) $(?:.F=.o)
	$(RANLIB) $(TDVRLIB)

sources: $(CPPSRCS)
	co -l RCS/*.F
