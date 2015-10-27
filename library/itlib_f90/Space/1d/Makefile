#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../../../itdvrlib.a
BLIULIB = ../../../bliulib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	h1onv.f \
	h12dvr.f \
	h1e.f \
	h12h0.f \
	h1v.f \
	rham1.f \
	setnd1.f \
	v12dvr.f \
	v12h0.f \
	v1e.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f=.o)

all: $(BLIULIB)
$(BLIULIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.f=.o)
