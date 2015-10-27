#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6055
#
FC = f77
FFLAGS = -c +OP

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6055/bfint1.f \
	../m6055/bfints.f \
	../m6055/bfprnt.f \
	../m6055/bndprn.f \
	../m6055/bsprnt.f \
	../m6055/crvadd.f \
	../m6055/ffbf.f \
	../m6055/ffint1.f \
	../m6055/ffints.f \
	../m6055/ffprnt.f \
	../m6055/filc12.f \
	../m6055/filcf1.f \
	../m6055/filcf2.f \
	../m6055/fildf1.f \
	../m6055/fildf2.f \
	../m6055/grdprn.f \
	../m6055/mkbes.f \
	../m6055/mkpow.f \
	../m6055/scalev.f \
	../m6055/tomobs.f \
	../m6055/vprnt.f \
	../m6055/yprnt.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6055

$(BINDIR)/m6055: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6055





