#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6287
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
	cdiis.f \
	chk.f \
	convt.f \
	diagx.f \
	diis.f \
	errvec.f \
	ham0.f \
	linsys.f \
	mkrhs.f \
	newply.f \
	prep.f \
	seidel.f \
	xmtrx.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6287

$(BINDIR)/m6287: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6287
