#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M5000
#
FC = fort

LD = fort
LIB2USE = orthlib
CPPFLAGS = -D$(LIB2USE)
FFLAGS = -c $(CPPFLAGS) -O5
LDFLAGS =
MAIN = /usr/home/bis/mesa
BINDIR = $(MAIN)/bin
MESALIB = $(MAIN)/library/mesalib.a
MDLIB= /usr/lib/libcxml.a
LIBDVR = $(MAIN)/library/PolyLib/$(LIB2USE).a
POTLIB = $(MAIN)/library/potlib/potlib.a


SRCS = \
	basfun.F \
	cmpre.F \
	chkfun.F \
	drvkmt.F \
	drvhmt.F \
	driver.F \
	fnlham.F \
	ovlp.F \
	renorm.F

.F.o:

	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m5000

$(BINDIR)/m5000: $(SRCS:.F=.o)	$(LIBDVR) $(POTLIB) $(MESALIB) $(MDLIB)
				$(LD) $(LDFLAGS) $(SRCS:.F=.o)  \
				$(LIBDVR) $(POTLIB) \
				$(MESALIB) $(MDLIB) -o $(BINDIR)/m5000
sources: $(SRCS)
	co -l RCS/*.F*
