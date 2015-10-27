#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6065
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6065/aoout.f \
	../m6065/fnlint.f \
	../m6065/intgl.f \
	../m6065/meopt.f \
	../m6065/moout.f \
	../m6065/newvec.f \
	../m6065/optpot.f \
	../m6065/ovrlp.f \
	../m6065/sclvlm.f \
	../m6065/tomo.f  \
	../m6065/zimag.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6065

$(BINDIR)/m6065: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6065





