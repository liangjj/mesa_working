#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6008
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6008/bbnum.f \
	../m6008/finbf.f \
	../m6008/finff.f \
	../m6008/kohnmt.f \
	../m6008/modden.f \
	../m6008/modnum.f \
	../m6008/newlst.f \
	../m6008/rdham.f \
	../m6008/rearrg.f \
	../m6008/resnum.f \
	../m6008/wrdeno.f \
	../m6008/wrnum.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6008

$(BINDIR)/m6008: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6008



