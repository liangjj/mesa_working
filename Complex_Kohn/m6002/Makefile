#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6002
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
	../m6002/c2v.f \
	../m6002/d2h.f \
	../m6002/dmatin.f \
	../m6002/fmch.f \
	../m6002/generf.f \
	../m6002/gfunct.f \
	../m6002/lgndre.f \
	../m6002/nlmxyz.f \
	../m6002/prntv.f \
	../m6002/rdbsis.f \
	../m6002/rdmom.f \
	../m6002/stuff0.f \
	../m6002/stuff1.f \
	../m6002/stuff2.f \
	../m6002/stuff3.f \
	../m6002/stuff4.f \
	../m6002/stuff5.f \
	../m6002/stuff6.f \
	../m6002/vecin.f \
	../m6002/vints.f \
	../m6002/vlr.f \
	../m6002/vstat.f \
	../m6002/wrtbas.f \
	../m6002/wrtpot.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6002

$(BINDIR)/m6002: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6002




