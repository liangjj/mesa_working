#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6236
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = /newdsk/mesa/bin
MESCTLIB = /newdsk/mesa/library/mesctlib.a
VECLIB = /usr/lib/libvec.a

GET = sccs get

SRCS = \
	egnpse.f \
	filpot.f \
	kmtrx.f \
	openc.f \
	potntl.f \
	rmprop.f \
	rmtrx.f \
	rmtxp.f \
	rprop1.f \
	transp.f \
	vcpl.f \
	wdcalc.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6236

$(BINDIR)/m6236: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6236
