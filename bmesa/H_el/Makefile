#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8020
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
LIBDIR = /usr1/mesa/library
BINDIR = /usr1/mesa/bin
MESCTLIB = $(LIBDIR)/mesctlib.a
VECLIB = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	drv.f \
	exnonit.f \
	fulsol.f \
	f1s.f \
	green.f \
	grid.f \
	hcoll.f \
	integ.f \
	local.f \
	nonit.f \
	phys.f \
	rdlam.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8020

$(BINDIR)/m8020: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m8020

