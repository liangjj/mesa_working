#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#

MAIN = /usr/home/bis/mesa
BINDIR = $(MAIN)/bin
MESALIB = $(MAIN)/library/mesalib.a
DVRLIB = $(MAIN)/library/dvrlib/dvrlib.a
ORTHLIB = $(MAIN)/library/orthlib/orthlib.a
POTLIB = $(MAIN)/library/potlib/potlib.a
MDLIB = /usr/lib/libcxml.a
FC = fort
LD = fort
LDFLAGS =
CPPFLAGS = -DDECPOINTER
FFLAGS = -c $(CPPFLAGS) -O4
RCS = RCS
CO = co -l
REL = 

CPPSRCS = \
	main.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m993
$(BINDIR)/m993: $(CPPSRCS:.F=.o) $(DVRLIB) $(ORTHLIB) $(POTLIB) \
				 $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(CPPSRCS:.F=.o) $(DVRLIB) $(ORTHLIB) \
					  $(POTLIB) $(MESALIB) $(MDLIB) \
					-o $(BINDIR)/m993

sources: $(CPPSRCS)
	$(CO) $(RCS)/*.F*

clean:
	rm *.o

