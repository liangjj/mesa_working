#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6003
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6003/fact.f \
	../m6003/legend.f \
	../m6003/miscfn.f \
	../m6003/prnty.f \
	../m6003/sph.f \
	../m6003/ylm.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6003

$(BINDIR)/m6003: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6003




