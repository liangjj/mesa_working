#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6004
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

SRCS = \
	../m6004/fitbes.f \
	../m6004/mkpow.f \
	../m6004/mkspln.f \
	../m6004/mkx.f \
	../m6004/prspln.f \
	../m6004/rcbesb.f \
	../m6004/rcbesf.f \
	../m6004/spline.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6004

$(BINDIR)/m6004: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6004





