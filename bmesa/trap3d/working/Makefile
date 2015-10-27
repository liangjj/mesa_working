#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7000
#
FC = f77
FFLAGS = -c 
LD = f77
LDFLAGS =
BINDIR = /usr1/mesa/bin
MESCTLIB = /usr1/mesa/library/mesctlib.a

GET = sccs get

SRCS = \
	addvnl.f \
	atmdat.f \
	chk.f \
	convt.f \
	diagx.f \
	diis.f \
	dropfn.f \
	dvd.f \
	dvr2h0.f \
	errmat.f \
	errfac.f \
	fock.f \
	fpsi.f \
	fxpt.f \
	hamil.f \
	ham0.f \
	h02dvr.f \
	honv.f \
	lschr.f \
	newerr.f \
	newfck.f \
	newfn.f \
	newvec.f \
	nlschr.f \
	nldvd.f \
	pkmat.f \
	potnl.f \
	prepfn.f \
	prtwfn.f \
	rdham.f \
	rdiag.f \
	setind.f \
	tcart.f \
	tcylin.f \
	tfsol.f \
	totpt.f \
	trap3d.f \
	trevod.f \
	trdiis.f \
	trmat.f \
	truth.f \
	tstovl.f \
	tthet.f \
	ua2d.f \
	ua3d.f \
	uta2d.f \
	uta3d.f \
	v2d.f \
	v3d.f \
	v2od.f \
	v3od.f \
	vmat.f \
	vmcosp.f \
	vpert.f \
	vtrial.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7000

$(BINDIR)/m7000: $(SRCS:.f=.o) $(MESCTLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) -o $(BINDIR)/m7000
