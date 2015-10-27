#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6296
#
FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS) 

LD = fort
LDFLAGS =

MAIN = /usr/home/bis/mesa
LIBDIR = $(MAIN)/library
BINDIR = $(MAIN)/bin
DVRLIB = $(LIBDIR)/dvrlib/dvrlib.a
COORDLIB = $(LIBDIR)/PolyLib/coordlib.a
POTLIB = $(LIBDIR)/potlib/potlib.a
HAMLIB = $(LIBDIR)/hamlib/hamlib.a
MESALIB = $(LIBDIR)/mesalib.a
PDELIB = $(LIBDIR)/pdelib/pdelib.a
DVDLIB = $(LIBDIR)/dvdson/dvdlib.a
ITDVRLIB  = $(LIBDIR)/itlib/itdvrlib.a
MDLIB = /usr/lib/libcxml.a


CPPSRCS = \
	dvddat.F \
	lindat.F \
	mgdvr.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6296

$(BINDIR)/m6296: $(CPPSRCS:F=o) $(DVRLIB) $(COORDLIB) $(HAMLIB) \
                               $(POTLIB) $(DVDLIB) $(ITDVRLIB) $(MESALIB) $(MDLIB)  
	$(LD) $(LDFLAGS) $(CPPSRCS:F=o) $(DVRLIB) $(COORDLIB) \
				       $(HAMLIB) $(POTLIB) $(DVDLIB) \
	                               $(ITDVRLIB) $(MESALIB) $(MDLIB) \
                                       -o $(BINDIR)/m6296
sources: $(CPPSRCS)
	co -l RCS/*.F*



