#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6018
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6018/bbkin.f \
	../m6018/filbf.f \
	../m6018/filff.f \
	../m6018/filmat.f \
	../m6018/finbf.f \
	../m6018/finff.f \
	../m6018/onemat.f \
	../m6018/ovfil.f \
	../m6018/tonwbf.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6018

$(BINDIR)/m6018: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6018



