#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6060
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6060/aoout.f \
	../m6060/ffbf.f \
	../m6060/filcf1.f \
	../m6060/fildf1.f \
	../m6060/mbes.f \
	../m6060/mkpow.f \
	../m6060/moout.f \
	../m6060/newlst.f \
	../m6060/optmat.f \
	../m6060/optpot.f \
	../m6060/ovrlp.f \
	../m6060/sclvlm.f \
	../m6060/tomo.f  \
	../m6060/zimag.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6060

$(BINDIR)/m6060: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6060





