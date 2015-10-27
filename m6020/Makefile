#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6020
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
	../m6020/drver.f \
	../m6020/fitbes.f \
	../m6020/fnlfun.f \
	../m6020/mkx.f \
	../m6020/mkxexp.f \
	../m6020/numerv.f \
	../m6020/potntl.f \
	../m6020/prbes.f \
	../m6020/prspln.f \
	../m6020/splinc.f \
	../m6020/splinr.f \
	../m6020/tolggr.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6020

$(BINDIR)/m6020: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6020




