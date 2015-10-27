#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m8001
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
	../onelam/conmat.f \
	../onelam/genvec.f \
	../onelam/grnset.f \
	../onelam/linslv.f \
	../onelam/m8001.f \
	../onelam/modpot.f \
	../onelam/onelin.f \
	../onelam/phase.f \
	../onelam/prrmat.f \
	../onelam/rhscal.f \
	../onelam/slvit.f \
	../onelam/vscale.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8001

$(BINDIR)/m8001: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m8001




