#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m1204
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../../../bin
MESCTLIB = ../../../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get

SRCS = \
	../old/bndryc.f \
	../old/bndryn.f \
	../old/colloc.f \
	../old/derivs.f \
	../old/filham.f \
	../old/funct.f \
        ../old/matpp.f \
	../old/matqp.f \
	../old/matqq.f \
        ../old/mkg.f \
	../old/mkgrd.f \
	../old/num3pt.f \
	../old/num5pt.f \
	../old/num7pt.f \
	../old/numerv.f \
	../old/phase.f \
        ../old/potntl.f \
	../old/veff.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1204

$(BINDIR)/m1204: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m1204




