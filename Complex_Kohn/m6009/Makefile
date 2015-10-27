#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6009
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6009/clvcmp.f \
	../m6009/cvbfr.f \
	../m6009/ffmat.f \
	../m6009/frmcmp.f \
	../m6009/kohnop.f \
	../m6009/matm.f \
	../m6009/mkfrbn.f \
	../m6009/mkfree.f \
	../m6009/rddeno.f \
	../m6009/rgslc.f \
	../m6009/slvlin.f \
	../m6009/tnonvr.f \
	../m6009/tvar.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6009

$(BINDIR)/m6009: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6009



