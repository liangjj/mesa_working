#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6200
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
	adgauss.f \
	class.f \
	fillf.f \
	gamfun.f \
	gaussq.f \
	gbslve.f \
	satshl.f \
	gbtql2.f \
	setsph.f \
	shells.f \
	trnsf.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6200

$(BINDIR)/m6200: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6200




