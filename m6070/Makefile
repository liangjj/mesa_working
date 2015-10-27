#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6070
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6070/mrgmt.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6070

$(BINDIR)/m6070: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6070





