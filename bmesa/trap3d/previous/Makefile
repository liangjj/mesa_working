#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7500
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = /newdsk/mesa/bin
MESCTLIB = /newdsk/mesa/library/mesctlib.a
VECLIB = /usr/lib/libvec.a

GET = sccs get

SRCS = \
	chk.f \
	convt.f \
	diagx.f \
	diis.f \
	ddiis.f \
	dnmat.f \
	dvd.f \
	dvr2h0.f \
	errmat.f \
	fock.f \
	hamil.f \
	ham0.f \
	h02dvr.f \
	honv.f \
	hnlonv.f \
	init.f \
	lschr1.f \
	lschr2.f \
	lschr3.f \
	newfn.f \
	newply.f \
	newvec.f \
	nlschr.f \
	nldvd.f \
	potnl.f \
	rdham.f \
	rdiag.f \
	setind.f \
	trap3d.f \
	trmat.f \
	tstovl.f \
	v2d.f \
	v3d.f \
	v2od.f \
	v3od.f \
	vmat.f \
	vmcosp.f \
	vtrial.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7500

$(BINDIR)/m7500: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m7500
