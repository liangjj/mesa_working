#   @(#)Makefile	1.4   8/2/91
#
# Makefile for tdvrlib
#
TDVRLIB = tdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O4
GET = sccs get
REL = 

SRCS = \
	addd.f \
	basis.f \
	cmpslv.f \
	eigslv.f \
	fildr.f \
	filvt.f \
	gpaket.f \
	hamfl1.f \
	hamfl2.f \
	hamfl3.f \
	hamt.f \
	hamtot.f \
	initpsi.f \
	lindat.f \
	mkh1d.f \
	mkh2d.f \
	mkh3d.f \
	mkpsi0.f \
	mkhpsi.f \
	pltwfn.f \
	setd2.f \
	setd3.f \
	setd4.f \
	set.f \
	vect0.f \
	zfill1.f \
	zfill2.f \
	zfill3.f \
	zproj.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(TDVRLIB)
$(TDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(TDVRLIB) $(?:.f=.o)
	$(RANLIB)
