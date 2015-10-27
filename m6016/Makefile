#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6016
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6016/tchank.f \
	../m6016/tout.f \
	../m6016/wrmat.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6016

$(BINDIR)/m6016: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6016



