#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6286
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
	errvec.f \
	fock.f \
	ham0.f \
	htree.f \
	init.f \
	linsys.f \
	newfn.f \
	mkrhs.f \
	newply.f \
	nlschr.f \
	prep.f \
	seidel.f \
	tdiis.f \
	tonew.f \
	tstcon.f \
	xmtrx.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6286

$(BINDIR)/m6286: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6286
