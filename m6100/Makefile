#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6100
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6100/basis.f \
	../m6100/fitbas.f \
	../m6100/freefn.f \
	../m6100/hbound.f \
	../m6100/hsmat.f \
	../m6100/kohn.f \
	../m6100/match.f \
	../m6100/modpot.f \
	../m6100/toorth.f \
	../m6100/tophse.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6100

$(BINDIR)/m6100: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6100





