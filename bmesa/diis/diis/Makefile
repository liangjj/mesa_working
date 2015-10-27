#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6288
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
	ham0.f \
	linsys.f \
	mkrhs.f \
	newply.f \
	precnd.f \
	prep.f \
	rdiag.f \
	rdiis.f \
	seidel.f \
	xmtrx.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6288

$(BINDIR)/m6288: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6288
