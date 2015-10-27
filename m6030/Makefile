#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6030
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

GET = sccs get

SRCS = \
	bbovlp.f \
	bbpot.f \
	bndbnd.f \
	onemat.f \
	tomobs.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6030

$(BINDIR)/m6030: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6030



