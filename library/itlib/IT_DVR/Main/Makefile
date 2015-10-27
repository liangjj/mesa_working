#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../../itdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

CPPFLAGS = -DDECPOINTER
#FFLAGS = -c -O5 $(CPPFLAGS)
FFLAGS = -c -O5 
REL = 

CPPSRCS = \
	dvderr.F \
	dvrdvd.F \
	hinit.F \
	htonv.F \
	hvdvr.F \
	hupdat.F \
	lindvr.F \
	linit.F \
	lres.F \
	lslv.F \
	newsol.F \
	prepit.F \
	set.F \
	trials.F \
	vonv.F \
	vtrblk.F \
	vtrsep.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(CPPSRCS:.F=.o) 
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.F=.o) 

sources: $(CPPSRCS)
	co -l RCS/*.F

