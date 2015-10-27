#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the pde routines
#
PDELIB = pdelib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	ab3d.f \
	atmdat.f \
	bfdat.f \
	cfill.f \
	chk.f \
	convt.f \
	dropfn.f \
	federin.f \
	filrpa.f \
	fromin.f \
	genrpa.f \
	gentrl.f \
	ham0.f \
	h0.f \
	h1.f \
	newfn.f \
	potnl.f \
	prepfn.f \
	rfill.f \
	rmtrx.f \
	tcart.f \
        tcylin.f \
        totpt.f \
        trmat.f \
	tthet.f \
        v2d.f \
        v2od.f \
        v3d.f \
        v3od.f \
	vdiag.f \
        vmat.f \
        vmat1d.f \
        vpert.f \
        wrdcnt.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(PDELIB)
$(PDELIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(PDELIB) $(?:.f=.o)



