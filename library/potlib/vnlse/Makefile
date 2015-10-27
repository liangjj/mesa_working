#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the potential library
#
POTLIB = ../potlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O4
GET = sccs get
REL = 

SRCS = \
	gfac.f \
	trgfac.f \
	vcprow.f \
	vecprd.f \
	vfac.f \
	vnlse.f \
	vtpert.f \
	vt0.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(POTLIB)
$(POTLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(POTLIB) $(?:.f=.o)

clean:
	rm *.o

