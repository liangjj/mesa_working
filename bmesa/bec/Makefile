#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6285
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
	cfcoef.f \
	chk.f \
	cmpslv.f \
	convt.f \
	diagx.f \
	ehamxt.f \
	fcoef.f \
	gschmt.f \
	ham0.f \
	ham1.f \
	ihamxt.f \
	mkpsi0.f \
	rdham.f \
	rdiag.f \
	rhs1.f \
	setind.f \
	tprop.f \
	tstsol.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6285

$(BINDIR)/m6285: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6285
