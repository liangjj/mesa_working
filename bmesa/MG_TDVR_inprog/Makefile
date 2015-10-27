#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8001
#
FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)

LD = fort
LDFLAGS =
LIBDIR = /usr/home/bis/mesa/library
BINDIR = /usr/home/bis/mesa/bin
DVRLIB = $(LIBDIR)/dvrlib/dvrlib.a
COORDLIB = $(LIBDIR)/PolyLib/coordlib.a
TDVRLIB = $(LIBDIR)/tdvrlib/tdvrlib.a
POTLIB = $(LIBDIR)/potlib/potlib.a
MESALIB = $(LIBDIR)/mesalib.a
PDELIB = $(LIBDIR)/pdelib/pdelib.a
DVDLIB = $(LIBDIR)/dvdson/dvdlib.a
ITDVRLIB = $(LIBDIR)/itlib/itdvrlib.a
MDLIB = /usr/lib/libcxml.a

GET = sccs get

SRCS = \
	psit.F \
	timdvr.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8001

$(BINDIR)/m8001: $(SRCS:.F=.o) $(TDVRLIB) $(DVRLIB) $(COORDLIB) $(POTLIB) $(ITDVRLIB) \
                               $(MESALIB) $(MDLIB)   
	$(LD) $(LDFLAGS) $(SRCS:.F=.o) $(TDVRLIB) $(DVRLIB) $(COORDLIB) $(POTLIB) \
                                $(ITDVRLIB) $(PDELIB) $(DVDLIB) \
				$(MESALIB) $(MDLIB) -o $(BINDIR)/m8001

