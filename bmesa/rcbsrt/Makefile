#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6234
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
BINDIR = /usr1/mesa/bin
MESCTLIB = /usr1/mesa/library/mesctlib.a
VECLIB = /usr1/mesa/library/libdxml.a

GET = sccs get

SRCS = \
	intrv.f \
	mkx.f \
	nwtrap.f \
	rcbsrt.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6234

$(BINDIR)/m6234: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6234
