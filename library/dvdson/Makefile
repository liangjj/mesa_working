#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
DVDLIB = dvdlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c 
GET = sccs get
REL = 

SRCS = \
	cdham.f \
	cdvd.f \
	cdvr2h0.f \
	ch02dvr.f \
	cindvd.f \
	clavec.f \
	clares.f \
	cnvec.f \
	creinit.f \
	crmres.f \
	ctrial.f \
	cua2d.f \
	cua3d.f \
	cuta2d.f \
	cuta3d.f \
	dvr2h0.f \
	eigord.f \
	frmres.f \
	h02dvr.f \
	hamil.f \
	honv.f \
	lares.f \
	lavec.f \
	lindvd.f \
	lsolve.f \
	newvec.f \
	rdham.f \
	rdiag.f \
	reinit.f \
	renorm.f \
	rpa.f \
	rpadiag.f \
	rpainit.f \
	rpaonv.f \
	rpavec.f \
	rsdvd.f \
	setind.f \
	tstovl.f \
	ua2d.f \
	ua3d.f \
	uta2d.f \
	uta3d.f \
	vtrial.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(DVDLIB)
$(DVDLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(DVDLIB) $(?:.f=.o)



