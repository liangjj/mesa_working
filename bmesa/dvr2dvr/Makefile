#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M4001
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
	drvlst.f \
	dvr2dvr.f \
	guesvc.f \
	ham.f \
	honv.f \
	iham.f \
	linslv.f \
	mkrhs.f \
	rdiag.f \
	rdham.f \
	setind.f \
	tr2ply.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m5001

$(BINDIR)/m5001: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m5001





























