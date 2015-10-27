#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6031
#
FC = f77
FFLAGS = -c 

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a

GET = sccs get

SRCS = \
	bsabox.f \
	diag.f \
	emake.f \
	kinmat.f \
	match.f \
	orthbs.f \
	poles.f \
	rmat.f \
	rmone.f \
	vmat.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6031

$(BINDIR)/m6031: $(SRCS:.f=.o)  ${MESCTLIB} $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  ${MESCTLIB} $(VECLIB) -o $(BINDIR)/m6031



