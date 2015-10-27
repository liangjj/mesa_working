#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6275
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
	davdag.f \
	david.f \
	diag.f \
	guess.f \
	honv.f \
	mkeham.f \
	tstdag.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6275

$(BINDIR)/m6275: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6275
