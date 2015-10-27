#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7102
#
FC = fort
FFLAGS = -c -O4

LD = fort
LDFLAGS =

MAIN = /usr/home/bis/mesa
BINDIR = $(MAIN)/bin
DVRLIB = $(MAIN)/library/dvrlib/dvrlib.a
POTLIB = $(MAIN)/library/potlib/potlib.a
HAMLIB = $(MAIN)/bmesa/R_Matrix/hamlib/hamlib.a
MESALIB = $(MAIN)/library/mesalib.a
PDELIB = $(MAIN)/library/pdelib/pdelib.a
DVDLIB = $(MAIN)/library/dvdson/dvdlib.a
MDLIB = /usr/lib/libcxml.a


GET = sccs get

SRCS = \
	conrmat.f \
	extrnl.f \
	hamio.f \
	kmtrx.f \
	mke.f \
	rmtrx2.f \
	relene.f \
	scat1d.f \
	setsct.f \
	smtrx.f \
	smtrx1.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7102

$(BINDIR)/m7102: $(SRCS:.f=.o) $(DVRLIB) $(HAMLIB) $(POTLIB) $(PDELIB) \
                               $(MESALIB) $(MDLIB)    
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(DVRLIB) $(HAMLIB) $(POTLIB) \
                                        $(MESALIB) $(PDELIB) $(MDLIB) \
					-o $(BINDIR)/m7102

