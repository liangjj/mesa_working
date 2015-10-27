#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6295
#
FC = fort
FFLAGS = -c -O4

LD = fort
LDFLAGS =
LIBDIR = /usr/home/bis/mesa/library
BINDIR = /usr/home/bis/mesa/bin
MESALIB = $(LIBDIR)/mesalib.a
PDELIB = $(LIBDIR)/pdelib/pdelib.a
DVDLIB = $(LIBDIR)/dvdson/dvdlib.a
ITDVRLIB = $(LIBDIR)/itlib/itdvrlib.a
MDLIB = /usr/lib/libcxml.a

GET = sccs get

SRCS = \
	addd.f \
	chkmem.f \
	cmpslv.f \
	dvderr.f \
	eigslv.f \
	fildr1.f \
	fildr2.f \
	fildr3.f \
	fillz.f \
	gett0.f \
	gi2gj.f \
	gpaket.f \
	hamfac.f \
	hamslv.f \
	hamfl1.f \
	hamfl2.f \
	hamfl3.f \
	hamt.f \
	hamtot.f \
	hamxyz.f \
	hupdat.f \
	idiag.f \
	i2j1d.f \
	i2j2d.f \
	i2j3d.f \
	i2j4d.f \
	initpsi.f \
	jacobi.f \
	lindat.f \
	mgdrv.f \
	mkh1d.f \
	mkh2d.f \
	mkh3d.f \
	mkpsi0.f \
	mkhpsi.f \
	newsol.f \
	nxtvec.f \
	pltwfn.f \
	pntmem.f \
	prmat.f \
	psit.f \
	qcklin.f \
	rdfun.f \
	rdke.f \
	rdpnt.f \
	setd2.f \
	setd3.f \
	setd4.f \
	setpnt.f \
	spamem.f \
	timdvr.f \
	vect0.f \
	vmat.f \
	vmat1d.f \
	vpert.f \
	vsmple.f \
	vt0.f \
	v1di2o.f \
	v2di2o.f \
	v3di2o.f \
	v4di2o.f \
	vivo.f \
	wrdadd.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6295

$(BINDIR)/m6295: $(SRCS:.f=.o) $(ITDVRLIB) $(PDELIB) $(DVDLIB) $(MESALIB) \
                               $(MDLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(ITDVRLIB) $(PDELIB) $(DVDLIB) \
                                        $(MESALIB)  $(MDLIB) -o $(BINDIR)/m6295

