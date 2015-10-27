#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M802
#
LD=$(FC)
LDFLAGS =

BINDIR = ../bin
MESALIB = $(MESA_LIB)/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m802/bkspac.f \
	../m802/bloksd.f \
	../m802/constd.f \
	../m802/drtcdd.f \
	../m802/drtind.f \
	../m802/extern.f \
	../m802/getchr.f \
	../m802/getcnt.f \
	../m802/getcod.f \
	../m802/getkey.f \
	../m802/getlin.f \
	../m802/getsym.f \
	../m802/intgrl.f \
	../m802/iprint.f \
	../m802/lwrwks.f \
	../m802/m802.f \
	../m802/movdrt.f \
	../m802/movref.f \
	../m802/mrgdrt.f \
	../m802/nextc.f \
	../m802/pakdrt.f \
	../m802/paldus.f \
	../m802/pandq.f \
	../m802/pm802.f \
	../m802/pqx.f \
	../m802/prtdrt.f \
	../m802/qntgrl.f \
	../m802/renumb.f \
	../m802/reordr.f \
	../m802/scan.f \
	../m802/skipbl.f \
	../m802/squash.f \
	../m802/srtdrt.f \
	../m802/tapesd.f \
	../m802/weight.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m802

$(BINDIR)/m802: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m802

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m802.f
	cat $(SRCS) > ../source/m802.f

print: link
	lpr ../source/m802.f

clean:
	rm -f *.o core



