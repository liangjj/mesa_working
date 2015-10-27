#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m8000
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
	../onelin/bndry.f \
	../onelin/ddrive.f \
	../onelin/genvec.f \
	../onelin/genvpt.f \
	../onelin/gen0.f \
	../onelin/linslv.f \
	../onelin/mkfx.f \
	../onelin/mkgrd.f \
	../onelin/mkrhs.f \
	../onelin/newtrm.f \
	../onelin/numerv.f \
	../onelin/onelin.f \
	../onelin/phase.f \
	../onelin/potntl.f \
	../onelin/slvit.f \
	../onelin/vdrive.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8000

$(BINDIR)/m8000: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m8000




