#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6090
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6090/basis.f \
	../m6090/conmat.f \
	../m6090/fitbas.f \
	../m6090/fitit.f \
	../m6090/freefn.f \
	../m6090/match.f \
	../m6090/modpot.f \
	../m6090/rhscal.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6090

$(BINDIR)/m6090: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6090





