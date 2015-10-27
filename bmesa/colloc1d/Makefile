#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6232
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../../bin
MESCTLIB = ../../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get

SRCS = \
	../colloc1d/bndryc.f \
	../colloc1d/bndryn.f \
	../colloc1d/colloc.f \
	../colloc1d/dtacol.f \
	../colloc1d/dtafd.f \
	../colloc1d/dtaft.f \
	../colloc1d/derivs.f \
	../colloc1d/filham.f \
	../colloc1d/funct.f \
        ../colloc1d/matpp.f \
	../colloc1d/matpq.f \
	../colloc1d/matqp.f \
	../colloc1d/matqq.f \
	../colloc1d/matslv.f \
        ../colloc1d/mkg.f \
	../colloc1d/mkgrd.f \
	../colloc1d/newqq.f \
	../colloc1d/num3pt.f \
	../colloc1d/num5pt.f \
	../colloc1d/num7pt.f \
	../colloc1d/numerv.f \
	../colloc1d/nwtprn.f \
	../colloc1d/phase.f \
        ../colloc1d/potntl.f \
        ../colloc1d/refine.f \
        ../colloc1d/solvep.f \
	../colloc1d/veff.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6232

$(BINDIR)/m6232: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6232




