#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6295
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
	cdiag.f \
	chk.f \
	clin.f \
	convt.f \
	diagx.f \
	dlin.f \
	drslv.f \
	dslv.f \
	ham0.f \
	h0slv.f \
	linslv.f \
	mkpsi0.f \
	preph.f \
	preslv.f \
	vmat.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6295

$(BINDIR)/m6295: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6295
