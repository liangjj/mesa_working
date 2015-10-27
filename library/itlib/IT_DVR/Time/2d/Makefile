#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../../../itdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	cham2.f \
	h2tonv.f \
	setnd3.f \
	setd3.f \
	v2t2dvr.f \
	v2t2h0.f \
	v2te.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f=.o)



