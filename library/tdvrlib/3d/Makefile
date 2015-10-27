#   @(#)Makefile	1.4   8/2/91
#
# Makefile for tdvrlib
#
TDVRLIB = ../tdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O4
GET = sccs get
REL = 

SRCS = \
	hamfl3.f \
	mkh3d.f \
	setd4.f \
	zfill3.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(TDVRLIB)
$(TDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(TDVRLIB) $(?:.f=.o)
	$(RANLIB) $(TDVRLIB)
