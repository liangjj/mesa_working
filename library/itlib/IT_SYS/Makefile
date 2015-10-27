#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
BLIULIB = ../bliulib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
CPPFLAGS = -DDECPOINTER

FFLAGS = -c -O5 $(CPPFLAGS)
GET = sccs get
REL = 

CPPSRCS = \
	bsmall.F \
	dvddat.F \
	honv.F \
	lares.F \
	lavec.F \
	lindvd.F \
	lsolve.F \
	rsdvd.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(BLIULIB)
$(BLIULIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.F=.o)

