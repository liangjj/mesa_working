#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M1200
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
	cmpare.f \
	intadd.f \
	mulgrd.f \
	relax.f \
	resid.f \
	resinj.f \
	rhoin.f \
	slvsml.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1200

$(BINDIR)/m1200: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m1200
