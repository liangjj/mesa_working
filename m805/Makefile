#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M805
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS = 

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../m805/bkspac.f \
	../m805/bloksd.f \
	../m805/constd.f \
	../m805/drtcdd.f \
	../m805/drtind.f \
	../m805/extern.f \
	../m805/getchr.f \
	../m805/getcnt.f \
	../m805/getcod.f \
	../m805/getkey.f \
	../m805/getlin.f \
	../m805/getsym.f \
	../m805/intgrl.f \
	../m805/iprint.f \
	../m805/lwrwks.f \
	../m805/m805.f \
	../m805/movdrt.f \
	../m805/movref.f \
	../m805/mrgdrt.f \
	../m805/nextc.f \
	../m805/pakdrt.f \
	../m805/paldus.f \
	../m805/pandq.f \
	../m805/pm805.f \
	../m805/pqx.f \
	../m805/prtdrt.f \
	../m805/qntgrl.f \
	../m805/renumb.f \
	../m805/reordr.f \
	../m805/scan.f \
	../m805/skipbl.f \
	../m805/squash.f \
	../m805/srtdrt.f \
	../m805/tapesd.f \
	../m805/weight.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m805
$(BINDIR)/m805: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m805

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m805.f
	cat $(SRCS) > ../source/m805.f

print: link
	lpr ../source/m805.f

clean:
	rm -f *.o core




