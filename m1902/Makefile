#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1902
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
MDLIB = $(MDLIB)
GET = sccs get
REL = 

SRCS = \
	../m1902/bf.f \
	../m1902/bfpntr.f \
	../m1902/del.f \
	../m1902/dirac.f \
	../m1902/efld.f \
	../m1902/eints.f \
	../m1902/fermi.f \
	../m1902/fmdel.f \
	../m1902/fmmult.f \
	../m1902/gcoord.f \
	../m1902/m1902.f \
	../m1902/mltpol.f \
	../m1902/pfour.f \
	../m1902/pone.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1902

$(BINDIR)/m1902: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1902

sources: $(SRCS)

clean:
	rm -f *.o core
$(SRCS):
	$(GET) $(REL) $@

