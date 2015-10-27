#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7010
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
	mkplot.f \
	plot.f \
	wfn.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7010

$(BINDIR)/m7010: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m7010
