#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6240
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
	cpoly.f \
	diff.f \
	fcoef.f \
	grid.f \
	gpoly.f \
	orpoly.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6240

$(BINDIR)/m6240: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6240
