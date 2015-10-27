#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m7000
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
	../m7000/cbasis.f \
	../m7000/drver.f \
	../m7000/fnlfun.f \
	../m7000/frelm.f \
	../m7000/genbes.f \
	../m7000/makefn.f \
	../m7000/mkbes.f \
	../m7000/mkx.f \
	../m7000/mkxexp.f \
	../m7000/nrmlz.f \
	../m7000/numerv.f \
	../m7000/potntl.f \
	../m7000/rbasis.f \
	../m7000/scalfn.f \
	../m7000/splinc.f \
	../m7000/splinr.f \
	../m7000/tolggr.f \
	../m7000/typylm.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7000

$(BINDIR)/m7000: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m7000




