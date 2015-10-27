#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6032
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
	diag.f \
	kinmat.f \
	tstnrm.f \
	vmat.f \
	wellp.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6032

$(BINDIR)/m6032: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6032



