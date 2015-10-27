#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6290
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
	buflin.f \
	chk.f \
	convt.f \
	diagx.f \
	dlin.f \
	drslv.f \
	dslv.f \
	ham0.f \
	honv.f \
	h0slv.f \
	ihamxt.f \
	lin.f \
	linslv.f \
	mkham.f \
	mkpsi0.f \
	newham.f \
	preph.f \
	preslv.f \
	rdiag.f \
	rdham.f \
	setind.f \
	vmat.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6290

$(BINDIR)/m6290: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6290
