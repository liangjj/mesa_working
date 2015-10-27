#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8010
#
FC = f90
FFLAGS = -c -O4

LD = f90
LDFLAGS =
MAIN = /usr1/mesa
LIBDIR = $(MAIN)/library
BINDIR = $(MAIN)/bin
DVRLIB = $(LIBDIR)/dvrlib/dvrlib.a
POTLIB = $(LIBDIR)/potlib/potlib.a
HAMLIB = $(LIBDIR)/hamlib/hamlib.a
MESCTLIB = $(LIBDIR)/mesctlib.a
PDELIB = $(LIBDIR)/pdelib/pdelib.a
DVDLIB = $(LIBDIR)/dvdson/dvdlib.a
ITDVRLIB  = $(LIBDIR)/itlib/itdvrlib.a
VECLIB = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	drvsplt.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8010

$(BINDIR)/m8010: $(SRCS:.f=.o)	$(DVRLIB) $(HAMLIB) $(POTLIB) $(DVDLIB) \
                                $(ITDVRLIB) $(PDELIB) $(MESCTLIB) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(DVRLIB) $(HAMLIB) $(POTLIB) \
                                        $(DVDLIB) $(ITDVRLIB) $(PDELIB) \
                                        $(MESCTLIB) $(VECLIB) \
					-o $(BINDIR)/m8010




