#  @(#)Makefile	1.8 8/1/91
#
# Makefile for M4550
#
FC = f77
FFLAGS = -c 

LD = f77
LDFLAGS =

BINDIR = /newdsk/mesa/bin
MESCTLIB = /newdsk/mesa/library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get
REL = 

SRCS = \
	bkwrd.f \
	cheby.f \
	forwrd.f \
	fulint.f \
	mkfun.f \
	mkgrd.f \
	pquad.f \
	quad.f \
	testnt.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m4550

$(BINDIR)/m4550: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m4550
