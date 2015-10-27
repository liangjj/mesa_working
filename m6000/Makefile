#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6000
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
	adgauss.f \
	gamfun.f \
	setsph.f \
	shells.f \
	trnsf.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6000

clean:
	rm -f *.o *~

$(BINDIR)/m6000: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6000



