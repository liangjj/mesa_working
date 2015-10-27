#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6001
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
	makecn.f \
	nlmxyz.f \
	orbs.f \
	rdbsis.f \
	wrtbas.f \
	wrtorb.f \
	wrttpe.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6001

$(BINDIR)/m6001: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6001



