#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6201
#
FC = f77
FFLAGS = -c +OP4
LD = f77
LDFLAGS =
BINDIR = /newdsk/mesa/bin
MESCTLIB = /newdsk/mesa/library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	decomp.f \
	fact.f \
	frms.f \
	legend.f \
	legint.f \
	miscfn.f \
	plmtof.f \
	scalfn.f \
	scmphi.f \
	sqrtwt.f \
	ylm.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6201

$(BINDIR)/m6201: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6201

