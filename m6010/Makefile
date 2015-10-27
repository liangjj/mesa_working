#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6010
#
FC = f77
FFLAGS = -c 
LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
GET = sccs get

SRCS = \
	basout.f \
	kohndt.f \
	newmat.f \
	trimc.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6010

$(BINDIR)/m6010: $(SRCS:.f=.o) $(MESCTLIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) -o $(BINDIR)/m6010


