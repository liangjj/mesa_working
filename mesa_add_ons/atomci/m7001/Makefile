#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m7001
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
	../m7001/bbden.f \
	../m7001/bcden.f \
	../m7001/bfden.f \
	../m7001/ccden.f \
	../m7001/cfden.f \
	../m7001/ffden.f \
	../m7001/lmtyp.f \
	../m7001/nprmax.f \
	../m7001/mkden.f \
	../m7001/mkfrs1.f \
	../m7001/mkonel.f \
	../m7001/mskone.f \
	../m7001/onebb.f \
	../m7001/onebc.f \
	../m7001/onebf.f \
	../m7001/onecc.f \
	../m7001/onecf.f \
	../m7001/oneff.f \
	../m7001/rhotst.f \
	../m7001/schmdt.f \
	../m7001/tomobb.f \
	../m7001/tomocc.f \
	../m7001/tomoff.f \
	../m7001/vearr.f \
	../m7001/vsarr.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7001

$(BINDIR)/m7001: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m7001




