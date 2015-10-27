#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M1202
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
	addint.f \
	interp.f \
	mulgrd.f \
	relax.f \
	resid.f \
	rhoin.f \
	rstrct.f \
	slvsml.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1202

$(BINDIR)/m1202: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m1202
